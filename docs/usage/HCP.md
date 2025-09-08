# High-Performance Computing (HPC)

## Did you know?

HippUnfold can automatically split its internal jobs and distribute them across a high-performance computing (HPC) cluster. This means you don’t have to rely only on your local workstation or manually batch subjects into separate job submissions. Instead, **Snakemake** handles job distribution, while HippUnfold provides sensible default resource settings for each step.

This setup lets you parallelize and speed up processing while still respecting resource allocation constraints defined by your cluster.

---

## Before you start

When running on a cluster, please check the following:

* **Installation path**: Make sure HippUnfold is installed in a directory that is accessible from all compute nodes (e.g., shared network storage).

* **Containers**: If you are using Apptainer/Singularity or Docker, confirm that the runtime is available on compute nodes.

* **Conda environments**: If using conda, ensure your environment directory is on a network-accessible path.

* **Cache directory**: By default, HippUnfold uses `~/.cache/hippunfold` (or your system’s default cache). If this is not visible to compute nodes, set a shared path, e.g.:

  ```bash
  export HIPPUNFOLD_CACHE_DIR=/YOUR/ACCESSIBLE/LOCATION/hippunfold/.cache
  ```

* **Internet access**: If compute nodes do not have internet access, pre-download all required data locally (e.g., by running HippUnfold on one subject locally first).

---

## SGE (Sun Grid Engine)

Sun Grid Engine (SGE) is a common cluster manager, typically called with commands such as:

```bash
qsub <your command>
```

with optional flags like `-q all.q` or `-q mica.q` (used in the [MICA Lab](https://mica-mni.github.io/)).

When running HippUnfold on SGE, replace the `--cores` option with `--profile` to point Snakemake to an SGE profile. For example, a profile might be stored at:

`~/.config/snakemake/sge/config.yml`

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: "qsub -cwd -q mica.q"
jobs: 20
cores: 16
max-threads: 16
latency-wait: 120
conda-frontend: mamba
conda-prefix: /data/mica3/.snakemake/conda
precommand: "source /host/cassio/export03/data/opt/miniconda3/etc/profile.d/conda.sh"
default-resources:
  - mem_mb=8000
  - time_h=6
```

⚠️ Adjust paths, memory/time defaults, and the `-q` option to suit your environment.

Submit HippUnfold with:

```bash
hippunfold <your usual arguments> --profile sge
```

---

## SLURM (Simple Linux Utility for Resource Management)

SLURM is another widely used cluster manager (e.g., on all **Compute Canada** systems). Jobs are usually submitted with:

```bash
sbatch <your command>
```

As with SGE, HippUnfold jobs can be distributed using a Snakemake profile for SLURM. A minimal example:

`~/.config/snakemake/slurm/config.yml`

```yaml
executor: slurm
jobs: 50
cores: 32
max-threads: 32
latency-wait: 120
conda-frontend: mamba
conda-prefix: /project/yourlab/.snakemake/conda
default-resources:
  - mem_mb=8000
  - time_h=6
cluster-cancel: scancel
cluster-status: sacct -j {jobid} --format=State --noheader
```

Submit HippUnfold with:

```bash
hippunfold <your usual arguments> --profile slurm
```

You may also need to customize `sbatch` options (e.g., partitions, QoS, GPUs). To do this, define a submission command in your profile, for example:

```yaml
cluster-generic-submit-cmd: "sbatch --partition=compute --time={resources.time_h}:00:00 --mem={resources.mem_mb}M"
```

---

## Summary

* Use `--profile` instead of `--cores` when submitting to a cluster.
* Configure your Snakemake profile (`sge/` or `slurm/`) to match your site’s environment.
* Ensure cache, conda environments, and containers are on network-accessible paths.
* Run one subject locally if needed to pre-download dependencies.

This setup will let HippUnfold scale efficiently across HPC systems with minimal manual job management.
