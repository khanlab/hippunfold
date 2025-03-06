FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="77000cfe577511b2d7412d0d6308c5f26bed08f08c2638d5ef1049a380cc4d9b"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/ants.yaml
#   prefix: /conda-envs/585d319787a6112fe84e87a9d19c0c71
#   name: ants
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - ants=2.5.4
#     - libitk=5.4.0
RUN mkdir -p /conda-envs/585d319787a6112fe84e87a9d19c0c71
COPY workflow/envs/ants.yaml /conda-envs/585d319787a6112fe84e87a9d19c0c71/environment.yaml

# Conda environment:
#   source: workflow/envs/c3d.yaml
#   prefix: /conda-envs/0737b7ab614d545aebcff2c27069c7c8
#   name: c3d
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - convert3d
#     - libitk=5.4.0
RUN mkdir -p /conda-envs/0737b7ab614d545aebcff2c27069c7c8
COPY workflow/envs/c3d.yaml /conda-envs/0737b7ab614d545aebcff2c27069c7c8/environment.yaml

# Conda environment:
#   source: workflow/envs/greedy.yaml
#   prefix: /conda-envs/caf2b760e4f7f4f36ef9717ea851a8e2
#   name: greedy
#   channels:
#     - khanlab
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - greedyreg
RUN mkdir -p /conda-envs/caf2b760e4f7f4f36ef9717ea851a8e2
COPY workflow/envs/greedy.yaml /conda-envs/caf2b760e4f7f4f36ef9717ea851a8e2/environment.yaml

# Conda environment:
#   source: workflow/envs/laynii.yaml
#   prefix: /conda-envs/d2c2bc1f56ea2545cd9d78799f019f63
#   name: laynii
#   channels:
#     - khanlab
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - laynii
RUN mkdir -p /conda-envs/d2c2bc1f56ea2545cd9d78799f019f63
COPY workflow/envs/laynii.yaml /conda-envs/d2c2bc1f56ea2545cd9d78799f019f63/environment.yaml

# Conda environment:
#   source: workflow/envs/neurovis.yaml
#   prefix: /conda-envs/bb76b1c086d16037586f4127e3d6474d
#   name: neurovis
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - nilearn>=0.8.1
#     - matplotlib=3.4.2
#     - pandas>=1.2.0,<1.3.0
#     - seaborn>=0.11.2
RUN mkdir -p /conda-envs/bb76b1c086d16037586f4127e3d6474d
COPY workflow/envs/neurovis.yaml /conda-envs/bb76b1c086d16037586f4127e3d6474d/environment.yaml

# Conda environment:
#   source: workflow/envs/niftyreg.yaml
#   prefix: /conda-envs/bbfc526e877e3476dc749a1ff3bf690b
#   name: niftyreg
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - niftyreg
RUN mkdir -p /conda-envs/bbfc526e877e3476dc749a1ff3bf690b
COPY workflow/envs/niftyreg.yaml /conda-envs/bbfc526e877e3476dc749a1ff3bf690b/environment.yaml

# Conda environment:
#   source: workflow/envs/nnunet.yaml
#   prefix: /conda-envs/d9280394e98f273aa715b9a733219e6f
#   name: nnunet
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - nnunet=1.7.1
RUN mkdir -p /conda-envs/d9280394e98f273aa715b9a733219e6f
COPY workflow/envs/nnunet.yaml /conda-envs/d9280394e98f273aa715b9a733219e6f/environment.yaml

# Conda environment:
#   source: workflow/envs/pyunfold.yaml
#   prefix: /conda-envs/ff23927388e3a494bfbbbd88d7d06a94
#   name: pyunfold
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - numpy=1.21.2
#     - astropy>=4.0,<4.3.1
#     - nibabel>=3.2.1
#     - scikit-fmm>=2022.03.26
#     - scipy=1.7.1
#     - pandas>=1.2.0,<1.3.0
RUN mkdir -p /conda-envs/ff23927388e3a494bfbbbd88d7d06a94
COPY workflow/envs/pyunfold.yaml /conda-envs/ff23927388e3a494bfbbbd88d7d06a94/environment.yaml

# Conda environment:
#   source: workflow/envs/pyvista.yaml
#   prefix: /conda-envs/839cf35b0d634cd136762cc2d1b9ab23
#   name: pyvista
#   channels:
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - numpy=1.21.2
#     - nibabel>=3.2.1
#     - pyvista>=0.44.2
#     - pygeodesic
#     - vtk
#     - scipy
RUN mkdir -p /conda-envs/839cf35b0d634cd136762cc2d1b9ab23
COPY workflow/envs/pyvista.yaml /conda-envs/839cf35b0d634cd136762cc2d1b9ab23/environment.yaml

# Conda environment:
#   source: workflow/envs/workbench.yaml
#   prefix: /conda-envs/6db957862240ca568f8d54608874aef6
#   name: workbench
#   channels:
#     - khanlab
#     - conda-forge
#     - defaults
#   dependencies:
#     - python=3.9
#     - dev-connectome-workbench
RUN mkdir -p /conda-envs/6db957862240ca568f8d54608874aef6
COPY workflow/envs/workbench.yaml /conda-envs/6db957862240ca568f8d54608874aef6/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/585d319787a6112fe84e87a9d19c0c71 --file /conda-envs/585d319787a6112fe84e87a9d19c0c71/environment.yaml && \
    conda env create --prefix /conda-envs/0737b7ab614d545aebcff2c27069c7c8 --file /conda-envs/0737b7ab614d545aebcff2c27069c7c8/environment.yaml && \
    conda env create --prefix /conda-envs/caf2b760e4f7f4f36ef9717ea851a8e2 --file /conda-envs/caf2b760e4f7f4f36ef9717ea851a8e2/environment.yaml && \
    conda env create --prefix /conda-envs/d2c2bc1f56ea2545cd9d78799f019f63 --file /conda-envs/d2c2bc1f56ea2545cd9d78799f019f63/environment.yaml && \
    conda env create --prefix /conda-envs/bb76b1c086d16037586f4127e3d6474d --file /conda-envs/bb76b1c086d16037586f4127e3d6474d/environment.yaml && \
    conda env create --prefix /conda-envs/bbfc526e877e3476dc749a1ff3bf690b --file /conda-envs/bbfc526e877e3476dc749a1ff3bf690b/environment.yaml && \
    conda env create --prefix /conda-envs/d9280394e98f273aa715b9a733219e6f --file /conda-envs/d9280394e98f273aa715b9a733219e6f/environment.yaml && \
    conda env create --prefix /conda-envs/ff23927388e3a494bfbbbd88d7d06a94 --file /conda-envs/ff23927388e3a494bfbbbd88d7d06a94/environment.yaml && \
    conda env create --prefix /conda-envs/839cf35b0d634cd136762cc2d1b9ab23 --file /conda-envs/839cf35b0d634cd136762cc2d1b9ab23/environment.yaml && \
    conda env create --prefix /conda-envs/6db957862240ca568f8d54608874aef6 --file /conda-envs/6db957862240ca568f8d54608874aef6/environment.yaml && \
    conda clean --all -y
