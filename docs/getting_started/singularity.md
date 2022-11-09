# Running HippUnfold with Singularity

## Pre-requisities:
 1. Singularity or Apptainer is installed on your system. For more info, see the detailed [apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
 2. The following command-line tools are installed:
  - wget
  - tar
 3. Sufficient disk-space in your `/tmp` folder (>30GB), as well as your working folder to store the container (~15GB), raw data, and processed data (~XXGB per subject) 
 4. Sufficient CPU and memory - the more you have, the faster it will run, but we recommend at least 8 CPU cores and 16GB memory.


### Downloading the test dataset

Download and extract a single-subject BIDS dataset for this test:

    wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
    tar -xvf hippunfold_test_data.tar

This will create a `ds002168/` folder with a single subject, that has a 
both T1w and T2w images:

```
ds002168/
├── dataset_description.json
├── README.md
└── sub-1425
    └── anat
        ├── sub-1425_T1w.json
        ├── sub-1425_T1w.nii.gz
        ├── sub-1425_T2w.json
        └── sub-1425_T2w.nii.gz

2 directories, 6 files
```

### Downloading HippUnfold

Next we will pull the container. The container can be pulled from dockerhub, but this requires a large amount of disk space in your /tmp folder, since it has to convert from a docker container to a singularity container. To avoid this, we provide a Dropbox link to the singularity container itself:

    wget https://www.dropbox.com/s/jtf6zyy0u8sc2k6/khanlab_hippunfold_v1.2.0.sif

If you really need to pull the container from docker hub, you can use the following command, but beware, it is more prone to errors and will take up lots of system resources (e.g. ~70GB of free disk space):

    singularity pull khanlab_hippunfold_v1.2.0.sif docker://khanlab/hippunfold:v1.2.0


Note: If you encounter any errors pulling the container from dockerhub, it may be because you are running 
out of disk space in your cache folders. Note, you can change these locations 
by setting environment variables, however, using a network file system for the folders may result in poor performance and/or errors e.g.:
    
    export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity


### Running HippUnfold

Run HippUnfold without any arguments to print the short help:

    singularity run -e khanlab_hippunfold_v1.2.0.sif 

Use the `-h` option to get a detailed help listing:

    singularity run -e khanlab_hippunfold_v1.2.0.sif -h

Note that all the Snakemake command-line options are also available in
HippUnfold, and can be listed with `--help-snakemake`:

    singularity run -e khanlab_hippunfold_v1.2.0.sif --help-snakemake


Now let's run it on the test dataset. The `--modality` flag is a 
required argument, and describes what image we use for segmentation. Here 
we will use the T1w image. We will also use the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.


    singularity run -e khanlab_hippunfold_v1.2.0.sif ds002168 ds002168_hippunfold participant -n --modality T1w


The first three arguments to HippUnfold (as with any BIDS App) are the input
folder, the output folder, and then the analysis level. The `participant` analysis 
level is used in HippUnfold for performing the segmentation, unfolding, and any
participant-level processing. The `group` analysis is used to combine subfield volumes
across subjects into a single tsv file.


When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    singularity run -e khanlab_hippunfold_v1.2.0.sif ds002168 ds002168_hippunfold participant -np --modality T1w | less


Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells HippUnfold how many cores to use.
 Using `--cores 8` means that HippUnfold will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (hippunfold on a single subject) may take ~30 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.


    singularity run -e khanlab_hippunfold_v1.2.0.sif ds002168 ds002168_hippunfold participant -p --cores all --modality T1w


Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/singularity_run.html) to ensure the container can read and write to yout input and output directories, respectively. You can bind paths easily by setting an 
environment variable, e.g. if you have a `/project` folder that contains your data, you can add it to the `SINGULARITY_BINDPATH` so it is available when you are running a container:

    export SINGULARITY_BINDPATH=/data:/data



After this completes, you should have a `ds002168_hippunfold` folder with outputs for the one subject.

If you alternatively want to run HippUnfold using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. In this case, since the T2w image in the 
test dataset has a limited FOV, we should also make use of the `--t1-reg-template` command-line option,
which will make use of the T1w image for template registration, since a limited FOV T2w template does not exist.

    singularity run -e khanlab_hippunfold_v1.2.0.sif ds002168 ds002168_hippunfold_t2w participant --modality T2w --t1-reg-template -p --cores all

Note that if you run with a different modality, you should use a separate output folder, since some of the files 
would be overwritten if not.




