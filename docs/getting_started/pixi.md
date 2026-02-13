# Running HippUnfold with pixi

HippUnfold can be installed and run using pixi on **Linux** system only. Pixi will manage all Python dependencies and non-python dependencies (c3d, greedy, ANTS) through conda environments.

**Note:** Pixi installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.


## For Users: Installing HippUnfold via Pixi


### Steps
 1. Install pixi (if not already installed):
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    ```
    
 2. Install hippunfold
    ```bash
    pixi global install hippunfold -c conda-forge -c khanlab -c bioconda

    ```


## Usage

### Test the installation

Run the following command to verify the installation:

```bash
hippunfold -h
```

You should see a help message listing all available command-line options.  
If this runs successfully, you’re ready to start processing data with HippUnfold!

## Running an example

You can try HippUnfold on a sample dataset to make sure everything works as expected.

First, download and extract a single-subject BIDS dataset for this test:

```bash
curl -L https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar -o hippunfold_test_data.tar
tar -xvf hippunfold_test_data.tar
```

This will create a `ds002168/` folder with a single subject that has both T1w and T2w images:

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

### Option 1: Run the full HippUnfold BIDS pipeline

Running HippUnfold using the T1w modality:

```bash
hippunfold ../ds002168 ../ds002168_hippunfold participant --modality T1w --cores all
```

This should run the full pipeline and place results in a new `ds002168_hippunfold/` folder.

### Option 2: Run HippUnfold for a single subject and modality directly

Alternatively, you can use the quick runner script to process just one image:

```bash
hippunfold-quick --input ../ds002168/sub-1425/anat/sub-1425_T1w.nii.gz --output ../ds002168_hippunfold_quick --modality T1w --subject 1425
```

This will run HippUnfold on the T1w image only and save outputs to the `ds002168_hippunfold_quick/` directory.

---

## Cache Directory

When running, HippUnfold automatically downloads and caches necessary resources such as atlases and templates to speed up subsequent runs.

By default, these are stored in the following directory:

```bash
~/.cache/hippunfold/
```

You can override this default cache location by setting the `HIPPUNFOLD_CACHE_DIR` environment variable:

```bash
export HIPPUNFOLD_CACHE_DIR=/path/to/custom/cache
```

This is useful when working on shared systems, when home directory storage is limited, or if you wish to isolate data per project or user.

    
## Development
For development work, use the development environment which includes additional tools like formatters and linters:

```bash
pixi install --environment dev
```

Quality checks can be run using:
```bash
pixi run --environment dev quality_check  # Check code formatting
pixi run  --environment dev quality_fix    # Fix code formatting
```   


---

Happy unfolding!
