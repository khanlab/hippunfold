# Command-line interface

## HippUnfold Command-line interface

 The following can also be seen by entering ``hippunfold -h`` into your terminal. 

These are all the required and optional arguments HippUnfold accepts in order to run flexibly on many different input data types and with many options, but in most cases only the required arguments are needed. 


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
module: snakemake
func: get_argument_parser
prog: snakemake
---
```


