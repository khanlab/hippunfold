# --- Stage 1: Build & Environment Prep ---
FROM ghcr.io/prefix-dev/pixi:latest AS build

WORKDIR /src
COPY . .

# Install pixi environment and create shell-hook
RUN pixi install --locked
RUN pixi shell-hook -s bash > /shell-hook
RUN mkdir -p /app && \
    echo "#!/bin/bash" > /app/entrypoint.sh && \
    cat /shell-hook >> /app/entrypoint.sh && \
    echo 'exec "$@"' >> /app/entrypoint.sh

# MOVE THE HEAVY GENERATION HERE (Inside Stage 1)
RUN set -e && \
    /app/entrypoint.sh hippunfold test_data/bids_singleT2w test_out participant --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_multiT2w test_out participant --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_T1w test_out participant --modality T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_hippb500 test_out participant --modality hippb500 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_singleT2w_longitudinal test_out participant --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_dsegtissue test_out participant --modality dsegtissue --derivatives test_data/bids_dsegtissue --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_singleT2w test_out participant --modality T2w --t1_reg_template --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_singleT2w test_out participant --modality T2w --output_space T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_singleT2w test_out participant --modality T2w --generate-myelin-map --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_T1w test_out participant --modality T1w --use-template-seg --template MBMv3 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    /app/entrypoint.sh hippunfold test_data/bids_dsegtissue test_out group_create_atlas --modality dsegtissue --derivatives test_data/bids_dsegtissue --new-atlas-name mytestatlas --new_atlas_subfields_from native --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    rm -rf /root/.cache /root/.conda/pkgs/* /src/.pixi/envs/default/pkgs/*

# --- Stage 2: Tiny Production Image ---
FROM ubuntu:24.04 AS production

WORKDIR /src

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy everything over from build stage (including pre-built conda-envs)
COPY --from=build /src /src
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh

ENV PYTHONNOUSERSITE=1
ENV SNAKEMAKE_PROFILE=/src/hippunfold/workflow/profiles/docker-conda

RUN echo '#!/bin/bash' > /usr/local/bin/hippunfold-quick && \
    echo 'exec /app/entrypoint.sh /src/hippunfold/run_quick.py "$@"' >> /usr/local/bin/hippunfold-quick && \
    chmod +x /usr/local/bin/hippunfold-quick

ENTRYPOINT ["/app/entrypoint.sh"]