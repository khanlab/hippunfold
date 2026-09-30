# Running on a SLURM cluster

```{note}
New in v2.2.0: HippUnfold includes resource specifications for every rule and groups short-running jobs together, so it can be run directly on a SLURM cluster with Snakemake's SLURM executor ([#574](https://github.com/khanlab/hippunfold/pull/574)).
```

HippUnfold is a [Snakemake](https://snakemake.readthedocs.io/) workflow, so it can submit its jobs to a SLURM scheduler instead of running them all on one machine. This is the recommended way to run HippUnfold on a cluster: you run `hippunfold` once for your whole dataset, and Snakemake submits each step as a SLURM job with appropriate memory, runtime and CPU requests, waits for them to finish, and submits the next steps as their inputs become available. Steps that need internet access (downloading the nnU-net model, templates and atlases) are run by the main `hippunfold` process itself rather than submitted, so they work even if your compute nodes have no internet access.

The [SLURM executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) (`snakemake-executor-plugin-slurm`) is installed with HippUnfold, so no additional installation is needed.[^executors]

## Quick start

```bash
hippunfold /path/to/bids /path/to/output participant --modality T1w \
    --executor slurm \
    --jobs 50 \
    --default-resources slurm_account=<your-account> \
    --retries 2
```

If your cluster does not require an account, replace `--default-resources slurm_account=<your-account>` with `--slurm-no-account`.

The `hippunfold` process must keep running until the workflow is complete, since it is what submits and monitors the jobs. Run it on a login node in a terminal multiplexer such as `tmux` or `screen`, or submit it as its own long, low-resource SLURM job (e.g. 1 CPU, 4 GB, for as long as you expect the whole dataset to take) if your cluster does not allow long-running processes on login nodes. In the latter case, make sure that job runs on a node with internet access, or that the models, templates and atlases have already been downloaded to the cache (see [](download-first)).

## Common options

| Option | Description |
|---|---|
| `--executor slurm` | Submit jobs to SLURM. |
| `--jobs N` | Maximum number of SLURM jobs submitted at the same time. |
| `--default-resources slurm_account=<acct>` | Account to charge jobs to. |
| `--slurm-no-account` | Submit without an account, for clusters that do not use them. |
| `--default-resources slurm_partition=<name>` | Partition to submit jobs to (otherwise the cluster's default partition is used). |
| `--retries N` | Resubmit failed jobs up to N times. Memory and runtime requests double with each attempt (see [below](#resources-and-retries)). |
| `--slurm-logdir <dir>` | Where to write the SLURM log files (default: `.snakemake/slurm_logs` in the output directory). |
| `--slurm-keep-successful-logs` | Keep SLURM logs of successful jobs, not just failed ones. |
| `--slurm-jobname-prefix <prefix>` | Prefix added to SLURM job names, e.g. `hippunfold`. |
| `--slurm-efficiency-report` | Write a report of CPU and memory efficiency of the jobs when the workflow finishes. |
| `--latency-wait N` | Seconds to wait for output files to appear on shared filesystems (increase if you see `MissingOutputException`). |

Several resources can be given to `--default-resources` at once, e.g. `--default-resources slurm_account=def-mylab slurm_partition=cpu`. The full list of SLURM options is in the [SLURM executor plugin documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html), or run `hippunfold --help-snakemake`.

### Using a profile

Instead of typing these options each time, you can put them in a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), e.g. `~/.config/snakemake/slurm/config.yaml`:

```yaml
executor: slurm
jobs: 50
retries: 2
latency-wait: 60
default-resources:
  slurm_account: def-mylab
  slurm_partition: cpu
```

and then run:

```bash
hippunfold /path/to/bids /path/to/output participant --modality T1w --profile slurm
```

## How jobs are submitted

Each HippUnfold rule specifies its memory (`mem_mb`), runtime and number of threads, based on profiling runs. Most steps take only seconds, so rather than submitting each one as its own SLURM job, consecutive short steps are combined into *group jobs*. For each subject, the workflow is submitted as roughly the following SLURM jobs:

| Group | Contents | Jobs per subject |
|---|---|---|
| `preproc` | Import, bias correction, registration to template, cropping to coronal-oblique space | 1 |
| `run_inference_*` | nnU-net segmentation (not grouped; the largest job, ~32 GB memory) | 1 per hemisphere |
| `qc_nnunet` | QC of the nnU-net segmentation | 1 per hemisphere |
| `shapeinject` | Template shape injection | 1 per hemisphere |
| `surf` | Laplace coordinates, surfaces, metrics, unfolded registration, subfields and QC | 1 per hemisphere and label (hipp, dentate) |
| `subj` | CIFTI and spec files combining both hemispheres, and subfield volumes | 1–2 |

For a typical T1w subject this is about 13 SLURM jobs. With the nnU-net inference on CPU, a subject takes roughly 15 minutes from start to finish when jobs start right away. Jobs for different subjects run in parallel, up to `--jobs` at a time.

The runtime requested for a group job is the sum of its steps' runtimes, so it is longer than the group actually takes (e.g. the `surf` groups request about an hour but usually finish in a few minutes).

(resources-and-retries)=
## Resources and retries

Memory and runtime requests scale with the attempt number: when a job is resubmitted with `--retries`, it requests double the memory and runtime of the previous attempt. Using `--retries 2` is recommended, so that jobs that run out of memory or time on unusually large images are automatically resubmitted with more.

You can also override the resources of individual rules, for example for high-resolution or ex vivo data:

```bash
hippunfold /path/to/bids /path/to/output participant --modality T2w \
    --executor slurm --jobs 50 --slurm-no-account \
    --set-resources laplace_beltrami:mem_mb=36000 laplace_beltrami:runtime=120 \
    --set-threads run_inference_nnunet_v1=8
```

Rule names can be listed with `hippunfold ... --list-rules`. For submitted jobs, the rule (or group) name is stored in the SLURM job comment, which you can show with e.g. `squeue --me -o "%.10i %.8T %.10M %k"`. To check how much of the requested resources your jobs actually used, run with `--slurm-efficiency-report`.

## Running many subjects

Run `hippunfold` once on the whole dataset: Snakemake submits the jobs for all subjects, up to `--jobs` at a time, and takes care of downloading the shared models, templates and atlases once.

If you instead run several `hippunfold` processes at the same time (e.g. one per `--participant-label` in a SLURM job array), they share the cache directory and can race each other downloading the same files. In that case, run the `download` analysis level once beforehand, as described in [](download-first).

[^executors]: Snakemake supports many other executors besides SLURM (e.g. other cluster schedulers and cloud platforms); see the [Snakemake plugin catalog](https://snakemake.github.io/snakemake-plugin-catalog/). These can be installed in the same environment as HippUnfold to enable them, and used with `--executor <name>`.
