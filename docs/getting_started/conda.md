# Running HippUnfold with Conda

HippUnfold can be installed and run using Conda on **Linux and macOS** systems.

**Note:** Conda installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.

---

## For Users: Installing HippUnfold via Conda

These steps are intended for **end users** who simply want to run HippUnfold for hippocampal segmentation and unfolding.

### 1. Install Conda (if not already installed)

Follow the instructions at the official Conda documentation:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

---

### 2. Create and activate a new Conda environment

```bash
conda create --name hippunfold-env {{ conda_channel }} -c conda-forge -c bioconda hippunfold
conda activate hippunfold-env
```

---

### 3. Test the installation

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
hippunfold ds002168 ds002168_hippunfold participant --modality T1w --cores all
```

This should run the full pipeline and place results in a new `ds002168_hippunfold/` folder.

### Option 2: Run HippUnfold for a single subject and modality directly

Alternatively, you can use the quick runner script to process just one image:

```bash
hippunfold-quick --input ds002168/sub-1425/anat/sub-1425_T1w.nii.gz --output ds002168_hippunfold_quick --modality T1w --subject 1425
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

## For Developers & Contributors

These steps are intended for people who want to contribute to the development of HippUnfold or explore its internals.

### 1. Clone the HippUnfold GitHub repository

```bash
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
```

---

### 2. Create and activate a new Conda environment

```bash
conda env create -f hippunfold/hippunfold-dev.yml
conda activate hippunfold-dev
```

---

### 3. Run the development version of HippUnfold

You can run HippUnfold directly from the source directory using:

```bash
./hippunfold/run.py -h
```

This should print out the available command-line arguments for the tool.  
You’re now set up for development and contribution!

---

## Troubleshooting

If you encounter issues while setting up HippUnfold via Conda:

- Make sure you’re using the latest Conda:
  ```bash
  conda update -n base -c defaults conda
  ```
- Double-check that your environment is activated (`conda activate hippunfold-env` or `hippunfold-dev`)
- Try creating a fresh environment if problems persist
- Search for similar issues or open a new one in the [GitHub issues](https://github.com/khanlab/hippunfold/issues) page

---

Happy unfolding!
