FROM condaforge/miniforge3:latest

WORKDIR /src/

# Copy your code
COPY . /src/

# Disable user site packages
ENV PYTHONNOUSERSITE=1

# Use bash for the following RUNs
SHELL ["/bin/bash", "-c"]

# ---- ONE SINGLE RUN ----
RUN set -e && \
    conda install -n base -c conda-forge mamba -y && \
    mamba create -y -n snakebids-env -c conda-forge -c bioconda snakebids unzip && \
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate snakebids-env && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_multiT2w test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_T1w test_out participant  --modality T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_hippb500 test_out participant  --modality hippb500 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_singleT2w_longitudinal test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_dsegtissue test_out participant  --modality dsegtissue --derivatives test_data/bids_dsegtissue --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --t1_reg_template --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --output_space T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_T1w test_out participant  --modality T1w --use-template-seg --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --generate-myelin-map --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    ./hippunfold/run.py test_data/bids_dsegtissue test_out group_create_atlas  --modality dsegtissue --derivatives test_data/bids_dsegtissue --new-atlas-name mytestatlas --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs --conda-frontend mamba && \
    conda clean --all -y && \
    rm -rf /opt/conda/pkgs /root/.cache && \
    cp /src/entrypoint_quick.sh /usr/local/bin/hippunfold-quick 

ENV SNAKEMAKE_PROFILE=/src/hippunfold/workflow/profiles/docker-conda

# Set entrypoint
ENTRYPOINT ["/src/entrypoint.sh"]


