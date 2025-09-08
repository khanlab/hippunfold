# Command-line interface


## HippUnfold Command-line interface

HippUnfold expects to parse all required files from all subjects in a BIDS directory. Thus a typical HippUnfold could simply be:

```bash
hippunfold /PATH/TO/YOUR/DATA /PATH/TO/OUTPUT participant --modality T1w
```

However, HippUnfold also has many flags that can be used to customize both the workflow and the input files being parsed. These can be seen by entering ``hippunfold -h`` into your terminal, which returns the following:

```{argparse}
---
filename: ../hippunfold/run.py
func: get_parser
prog: hippunfold
---
```


## Snakemake command-line interface

In addition to the above command-line arguments, Snakemake arguments are also be passed at the `hippunfold` command-line. 

The most critical of these is the `--cores` or `-c` argument, which is a **required** argument for HippUnfold. 

The complete list of [Snakemake](https://snakemake.readthedocs.io/en/stable/) arguments are below, and mostly act to determine your environment and App behaviours. They will likely only need to be used for running in cloud environments or troubleshooting. These can be listed from the command-line with `hippunfold --help-snakemake`.  

```{argparse}
---
module: snakemake.cli
func: get_argument_parser
prog: snakemake
---
```
