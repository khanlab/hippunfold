# Running HippUnfold with Conda

HippUnfold can be installed and run using Conda on **Linux and macOS** systems.

**Note:** Conda installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.

---

## For Users: Installing HippUnfold via Conda

These steps are intended for **end users** who simply want to run HippUnfold for hippocampal segmentation and unfolding.

### 1. Install Miniconda (if not already installed)

Follow the instructions at the official Miniconda site:  
🔗 [https://www.anaconda.com/docs/getting-started/miniconda/install](https://www.anaconda.com/docs/getting-started/miniconda/install)

---

### 2. Create and activate a new Conda environment

```bash
conda create -n hippunfold-env python=3.10
conda activate hippunfold-env
```

---

### 3. Install HippUnfold from Bioconda

```bash
conda install -c conda-forge -c bioconda hippunfold
```

---

### 4. Test the installation

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
wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
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

Now let’s run HippUnfold using the T1w modality:

```bash
hippunfold ds002168 ds002168_hippunfold participant --modality T1w --cores all
```

This should run the full pipeline and place results in a new `ds002168_hippunfold/` folder.

---

## For Developers & Contributors

These steps are intended for people who want to contribute to the development of HippUnfold or explore its internals.

### 1. Create and activate a new Conda environment

```bash
conda create -n hippunfold-dev python=3.10
conda activate hippunfold-dev
```

---

### 2. Install Snakebids and development dependencies

```bash
conda install -c conda-forge -c bioconda snakebids
```

---

### 3. Clone the HippUnfold GitHub repository

```bash
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
```

---

### 4. Run the development version of HippUnfold

You can run HippUnfold directly from the source directory using:

```bash
python hippunfold/run.py -h
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
For more information, visit the [official documentation]( hippunfold.readthedocs.io ).
