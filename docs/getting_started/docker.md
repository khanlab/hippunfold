# Running HippUnfold with Docker on Windows

Note: These instructions assume you have Docker installed already on a Windows system. Docker can also run on Linux or MacOS with similar commands, but here we will assume the default Windows CLI is being used.

## First time setup

Open your Windows Command Prompt by clicking the bottom left `Windows` button and type `cmd` followed by `Enter`. This is where you will enter your HippUnfold commands. Feel free to make a new directory with `mkdir` or move to a directory you would like to work out of with `cd`, and for this example we will work from:

    cd c:\Users\jordan\Downloads\

Pull the container (this will take some time and storage space, but like an installation it only needs to be done once and can then be run on many datasets):

    docker pull khanlab/hippunfold:latest

Run HippUnfold without any arguments to print the short help:

    docker run -it --rm khanlab/hippunfold:latest    

Use the `-h` option to get a detailed help listing:

    docker run -it --rm khanlab/hippunfold:latest -h

Note that all the Snakemake command-line options are also available in
HippUnfold, and can be listed with `--help-snakemake`:

    docker run -it --rm khanlab/hippunfold:latest --help-snakemake

## Running an example

Download and extract a single-subject BIDS dataset for this test from [hippunfold_test_data.tar](https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar). Here we will also assume you chose to save and extract to the directory `c:\Users\jordan\Downloads\`.

This contains a `ds002168/` directory with a single subject, that has a both T1w and T2w images. 

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

Now let's run HippUnfold on the test dataset. Docker will need read/write access to the input and output directories, respectively. This is achieved with the `-v` flag. This 'binds' or 'mounts' a directory to a new directory inside the container.

    docker run -it --rm -v c:\Users\jordan\Downloads\ds002168:/bids -v c:\Users\jordan\Downloads\ds002168_hippunfold:/output khanlab/hippunfold:latest /bids /output participant --modality T1w -n

Explanation: 

`-v c:\Users\jordan\Downloads\ds002168:/bids` tells Docker to mount the directory `c:\Users\jordan\Downloads\ds002168` into a new directory inside the container named `/bids`. We then do the same things for our output directory named `ds002168_hippunfold`, which we mount to `/output` inside the container. These arguments are not specific to HippUnfold but rather are general ways to use Docker. You may want to familiarize yourself with [Docker options](https://docs.docker.com/engine/reference/run/).

Everything after we specified the container (`khanlab/hippunfold:latest`) are arguments to HippUnfold itself. The first of these arguments (as with any BIDS App) are the input directory (`/bids`), the output directory (`/output`), and then the analysis level (`participant`). The `participant` analysis 
level is used in HippUnfold for performing the segmentation, unfolding, and any
participant-level processing. The `group` analysis is used to combine subfield volumes
across subjects into a single tsv file.  The `--modality` flag is also required, and describes which image we use for segmentation. Here we used the T1w image. We also used the `--dry-run/-n`  option to just print out what would run, without actually running anything.

When you run the above command, a long listing will print out, describing all the rules that 
will be run. Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells HippUnfold how many cores to use.
 Using `--cores 8` means that HippUnfold will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (hippunfold on a single subject) may take ~30 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.

    docker run -it --rm -v c:\Users\jordan\Downloads\ds002168:/bids -v c:\Users\jordan\Downloads\ds002168_hippunfold:/output khanlab/hippunfold:latest /bids /output participant --modality T1w -p --cores all

After this completes, you should have a `ds002168_hippunfold` directory with outputs for the one subject.

## Exploring different options

If you alternatively want to run HippUnfold using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. In this case, since the T2w image in the 
test dataset has a limited FOV, we should also make use of the `--t1-reg-template` command-line option,
which will make use of the T1w image for template registration, since a limited FOV T2w template does not exist.

    docker run -it --rm -v c:\Users\jordan\Downloads\ds002168:/bids -v c:\Users\jordan\Downloads\ds002168_hippunfold_t2w:/output khanlab/hippunfold:latest /bids /output participant --modality T2w --t1-reg-template -p --cores all

Note that if you run with a different modality, you should use a separate output directory, since some of the files 
would be overwritten if not.



