Command line interface
--------------------
The following can also be seen by entering ``hippunfold -h`` into your terminal. 

These are all the required and optional arguments HippUnfold accepts in order to run flexibly on many different input data types and with many options, but in most cases only the required arguments are needed. 

For further Snakemake options including environment variables, cloud computing options, or other behaviours see `snakemake arguments <https://hippunfold.readthedocs.io/en/latest/usage/snakemake_cli.html>`_ (or ``hippunfold --help_snakemake`` in terminal), which are provided after HippUnfold arguments. 

.. argparse::
   :filename: ../hippunfold/run.py
   :func: get_parser
   :prog: hippunfold


