# Running HippUnfold with Docker

Note: These instructions assume you have Docker installed already on your system.


Download and extract a single-subject BIDS dataset for this test:

    wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
    tar -xvf hippunfold_test_data.tar

This will create a `ds002168/` folder with a single subject, that has a 
both T1w and T2w images. 

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


Pull the container:

    docker pull khanlab/hippunfold:latest

Run HippUnfold without any arguments to print the short help:

    docker run -it --rm \
    khanlab/hippunfold:latest    

Use the `-h` option to get a detailed help listing:

    docker run -it --rm \
    khanlab/hippunfold:latest \
    -h

Note that all the Snakemake command-line options are also available in
HippUnfold, and can be listed with `--help-snakemake`:

    docker run -it --rm \
    khanlab/hippunfold:latest \
    --help-snakemake


Now let's run it on the test dataset. The `--modality` flag is a 
required argument, and describes what image we use for segmentation. Here 
we will use the T1w image. We will also use the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -n


For those not familiar with Docker, the first three lines of this
example are generic Docker arguments to ensure it is run with the safest
options and has permission to access your input and output directories
The fourth line specifies the HippUnfold Docker container, and the fifth line contains the required
arguments for HippUnfold, after which you can additionally specify optional arguments. You may want to familiarize yourself with
[Docker options](https://docs.docker.com/engine/reference/run/), and an
overview of HippUnfold arguments is provided in the [Command line
interface](https://hippunfold.readthedocs.io/en/latest/usage/cli.html)
documentation section.


The first three arguments to HippUnfold (as with any BIDS App) are the input
folder, the output folder, and then the analysis level. The `participant` analysis 
level is used in HippUnfold for performing the segmentation, unfolding, and any
participant-level processing. The `group` analysis is used to combine subfield volumes
across subjects into a single tsv file.


When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -np | less


Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells HippUnfold how many cores to use.
 Using `--cores 8` means that HippUnfold will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (hippunfold on a single subject) may take ~30 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T1w -p --cores all


After this completes, you should have a `ds002168_hippunfold` folder with outputs for the one subject.

If you alternatively want to run HippUnfold using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. In this case, since the T2w image in the 
test dataset has a limited FOV, we should also make use of the `--t1-reg-template` command-line option,
which will make use of the T1w image for template registration, since a limited FOV T2w template does not exist.

    docker run -it --rm \
    -v ds002168:/bids:ro \
    -v ds002168_hippunfold_t2w:/output \
    khanlab/hippunfold:latest \
    /bids /output participant --modality T2w --t1-reg-template -p --cores all

Note that if you run with a different modality, you should use a separate output folder, since some of the files 
would be overwritten if not.



