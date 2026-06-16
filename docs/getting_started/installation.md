# Installation (advanced)


These are more advanced ways to run hippunfold. To simply get started (install using pixi instead).

* **Conda (Linux/macOS)** — more flexible, good for HPC environments, developers.
* **Containers (Docker or Apptainer/Singularity)** — easiest, best for Windows, or when you want fully pinned dependencies.
  > **Note:** Containers must be given access to input/output directories by mounting them.

| Method                  | OS                            | Edit code? | Cluster profiles | Notes                                   |
| ----------------------- | ----------------------------- | ---------- | ---------------- | --------------------------------------- |
| Conda                   | Linux, macOS                  | ✅          | ✅                | Simple, container‑free. Not on Windows. |
| Docker                  | Windows, Linux, macOS (Intel) | ❌          | ❌                | Mount paths required; portable.         |
| Apptainer (Singularity) | Linux                         | ❌          | ❌                | Common on HPC; portable single `.sif`.  |


## Quickstart (choose one)

### Option A — Conda (Linux/macOS)

```bash
# create env & install (channels templated)
conda create --name hippunfold-env {{ conda_channel }} -c conda-forge -c bioconda hippunfold
conda activate hippunfold-env

# check it works
hippunfold -h
```

Run a one‑subject example:

```bash
curl -L https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar -o hippunfold_test_data.tar
tar -xvf hippunfold_test_data.tar
hippunfold ds002168 ds002168_hippunfold participant --modality T1w --cores all
```

### Option B — Containers (Docker or Apptainer)

Containers share nearly identical commands. Replace `docker run` with `apptainer run` depending on your environment. Both need directory mounts to access your data.

#### 1. Pull the container

Docker:

```bash
docker pull khanlab/hippunfold:{{ current_tag }}
```

Apptainer:

```bash
apptainer pull hippunfold_{{ current_tag }}.sif docker://khanlab/hippunfold:{{ current_tag }}
```

#### 2. Dry‑run

```bash
# Docker
docker run -it --rm \
  -v /directory/to/mount:/data \
  khanlab/hippunfold:{{ current_tag }} \
  /data/ds002168 /data/ds002168_hippunfold participant --modality T1w -n

# Apptainer
apptainer run -e \
  --bind /directory/to/mount:/data \
  hippunfold_{{ current_tag }}.sif \
  /data/ds002168 /data/ds002168_hippunfold participant --modality T1w -n
```

#### 3. Run with all cores

```bash
# Docker
docker run -it --rm \
  -v /directory/to/mount:/data \
  khanlab/hippunfold:{{ current_tag }} \
  /data/ds002168 /data/ds002168_hippunfold participant --modality T1w -p --cores all

# Apptainer
apptainer run -e \
  --bind /directory/to/mount:/data \
  hippunfold_{{ current_tag }}.sif \
  /data/ds002168 /data/ds002168_hippunfold participant --modality T1w -p --cores all
```

#### Note on mounting

* Docker uses `-v host_path:container_path`
* Apptainer uses `--bind host_path:container_path`

To avoid repeating mounts on every command, you can set:

```bash
# Apptainer
env APPTAINER_BINDPATH=/project:/project,/data:/data

# Docker (in docker-compose or wrapper scripts)
# e.g., export DOCKER_OPTS="-v /project:/project -v /data:/data"
```

This way, input/output paths are automatically accessible.

## Models & cache (Note)

As of **v1.3.0+**, containers no longer ship all models; they’re downloaded on demand.

Default cache:

```bash
~/.cache/hippunfold
```

Override with:

```bash
export HIPPUNFOLD_CACHE_DIR=/path/to/cache
```
