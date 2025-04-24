# Use a minimal Conda base image
FROM continuumio/miniconda3:latest

WORKDIR /src/

COPY . /src/

# avoid pre-downloading the models to make for lighter container
# ENV HIPPUNFOLD_CACHE_DIR=/opt/hippunfold_cache
ENV PYTHONNOUSERSITE=1

# (Optional) Update Conda and install Mamba for faster env creation
RUN conda update -n base -c conda-forge conda -y && \
    conda install -n base -c conda-forge mamba -y

# Install snakebids in a bootstrap Conda env (for running the Snakefile)
RUN mamba create -n snakebids-env -c conda-forge -c bioconda snakebids unzip -y && \
    echo "conda activate snakebids-env" >> ~/.bashrc


# Activate snakebids env and pre-create rule-specific envs
RUN bash -c "source ~/.bashrc && \
    conda activate snakebids-env && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --hemi R --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --hemi L --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_multiT2w test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_T1w test_out participant  --modality T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_hippb500 test_out participant  --modality hippb500 --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-env && \
    ./hippunfold/run.py test_data/bids_T1w_longitudinal test_out participant  --modality T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w_longitudinal test_out participant  --modality T2w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_dsegtissue test_out participant  --modality dsegtissue --derivatives test_data/bids_dsegtissue --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_dsegtissue_1hemi test_out participant  --modality dsegtissue --derivatives test_data/bids_dsegtissue_1hemi --hemi L --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --t1_reg_template --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --output_space T1w --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_T1w test_out participant  --modality T1w --use-template-seg --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_singleT2w test_out participant  --modality T2w --generate-myelin-map --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs && \
    ./hippunfold/run.py test_data/bids_dsegtissue test_out group_create_atlas  --modality dsegtissue --derivatives test_data/bids_dsegtissue --new-atlas-name mytestatlas --use-conda --conda-create-envs-only --cores all --conda-prefix /src/conda-envs"

# Run the pipeline
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "snakebids-env", "./hippunfold/run.py"]
CMD ["--use-conda", "--conda-prefix /src/conda-envs"]