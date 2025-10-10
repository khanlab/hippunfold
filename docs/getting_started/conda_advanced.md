# Conda (Advanced)

HippUnfold via Conda on **Linux/macOS**.

> Windows users: use containers instead (see Quickstart page).

These steps are intended for **developers or contributors** modifying HippUnfold code, or using advanced Snakemake features (e.g. cluster execution profiles).

### 1. Clone the repository

```bash
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
```

### 2. Create and activate a dev environment

Option A — Conda:

```bash
conda env create -f hippunfold/hippunfold-dev.yml
conda activate hippunfold-dev
```

Option B — Poetry (preferred for development):

```bash
# install poetry if not present (see https://python-poetry.org/docs/)
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
poetry install

# run via poetry
poetry run hippunfold
# or enter a shell
poetry shell
hippunfold
```

### 3. Code quality & formatting

We use [`poethepoet`](https://github.com/nat-n/poethepoet) with `black` and `snakefmt`.

```bash
poetry run poe quality_check
poetry run poe quality_fix
```

Inside a poetry shell, drop `poetry run`.

### 4. Dry‑run workflow tests

Use Snakemake’s dry‑run (`-n`) to verify changes. Example with fake test data included in the repo:

```bash
hippunfold test_data/bids_singleT2w test_out participant --modality T2w --use-singularity -np
```

Or run Snakemake directly:

```bash
cd hippunfold
snakemake -np
```

### 5. Wet‑run testing

Opensource test datasets are available on OSF: [https://osf.io/k2nme/](https://osf.io/k2nme/)

These cover different modalities, resolutions, and species. Use them to validate workflow changes (especially templates or modality‑specific logic).

### 6. Compute Canada / CBS clusters (optional)

For contributors working on HPC clusters (e.g. Graham, CBS):

* Use `pip` or `poetry` installs.
* Configure [neuroglia-helpers](https://github.com/khanlab/neuroglia-helpers) and [cc-slurm](https://github.com/khanlab/cc-slurm) profiles.
* Set env vars for shared model/cache directories:

  ```bash
  export HIPPUNFOLD_CACHE_DIR=/project/.../hippunfold_models
  export SNAKEMAKE_SINGULARITY_DIR=/project/.../containers
  ```
* Always run jobs inside `/localscratch` or `$SCRATCH`.

### 7. Models & cache

* Models are **downloaded automatically** as of v1.3.0+.
* Default cache: `~/.cache/hippunfold` (override with `HIPPUNFOLD_CACHE_DIR`).
* To pre‑download models (for offline compute nodes):

  ```bash
  hippunfold BIDS_DIR OUT_DIR participant --modality T1w --until download_model -c 1
  ```

### 8. Singularity/Apptainer cache overrides

```bash
export APPTAINER_CACHEDIR=/path/to/.cache/apptainer
export HIPPUNFOLD_CACHE_DIR=/path/to/.cache/hippunfold
```

When running with `--use-singularity`, containers are downloaded into `.snakemake` by default (override with `--singularity-prefix`).

---

## Troubleshooting

* Update Conda:

  ```bash
  conda update -n base -c defaults conda
  ```
* Ensure your env is active.
* Try a fresh env if problems persist.
* Check GitHub issues: [https://github.com/khanlab/hippunfold/issues](https://github.com/khanlab/hippunfold/issues)
