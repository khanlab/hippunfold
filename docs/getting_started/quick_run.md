# Quick Run Mode

HippUnfold includes a convenient `hippunfold_quick` command for processing single image files with improved efficiency through the use of temporary directories.

## Overview

The `hippunfold_quick` script is designed to:
- Accept single NIfTI image files (no BIDS structure required)
- Automatically create a temporary BIDS dataset
- Run the workflow in temporary directories (beneficial for network storage scenarios)
- Copy results back to your desired output location

This is particularly useful when:
- Your output directory is on a network drive (NFS, CIFS, etc.) where local temporary storage offers better I/O performance
- You want to quickly process a single image without manually creating a BIDS directory structure
- You need isolation between the workflow execution and final output location

## Installation

The `hippunfold_quick` command is installed automatically with HippUnfold:

```bash
pip install hippunfold
```

Or for development:

```bash
git clone https://github.com/khanlab/hippunfold
cd hippunfold
pip install -e .
```

## Basic Usage

Process a single T2w image:

```bash
hippunfold_quick sub-001_T2w.nii.gz /path/to/output
```

This will:
1. Create a temporary BIDS dataset with subject ID "001" (default)
2. Run hippunfold in temporary directories
3. Copy the results to `/path/to/output/hippunfold/sub-001/`
4. Clean up temporary files

## Options

### Subject ID

Specify a custom subject identifier:

```bash
hippunfold_quick brain.nii.gz /output --subject 123
```

### Modality

Specify the image modality (default: T2w):

```bash
hippunfold_quick brain.nii.gz /output --modality T1w
```

Supported modalities:
- `T1w` - T1-weighted image
- `T2w` - T2-weighted image
- `hippb500` - Hippocampal B500 DWI image

### Keep Temporary Files

Keep temporary directories for debugging:

```bash
hippunfold_quick input.nii.gz /output --keep-temp
```

### Passing Options to HippUnfold

All additional arguments are passed directly to hippunfold:

```bash
# Dry run with specific hemisphere
hippunfold_quick input.nii.gz /output -n --cores 8 --hemi L

# Use specific atlas and output density
hippunfold_quick input.nii.gz /output --atlas multihist7 --output_density 1mm
```

## How It Works

### Workflow

```
Input File → Temp BIDS Dataset → Temp Output Dir → Final Output Dir
   (you)       (auto-created)     (auto-created)    (copied results)
```

### Temporary Directory Structure

The script creates a temporary directory with the following structure:

```
/tmp/hippunfold_quick_XXXXX/
├── bids_input/
│   ├── dataset_description.json
│   └── sub-{subject}/
│       └── anat/
│           ├── sub-{subject}_{modality}.nii.gz
│           └── sub-{subject}_{modality}.json
└── hippunfold_output/
    ├── hippunfold/
    │   └── sub-{subject}/  ← This gets copied to final output
    ├── logs/
    └── work/
```

### Performance Benefits

When output directories are on network storage:
- Faster I/O for intermediate workflow files
- Reduced network traffic during processing
- Only final results are transferred over the network
- More reliable handling of temporary files

## Examples

### Process T2w image with dry run

```bash
hippunfold_quick ~/data/hippocampus_T2w.nii.gz /mnt/network/output -n
```

### Process T1w image with 16 cores

```bash
hippunfold_quick scan.nii.gz /output --modality T1w --cores 16
```

### Process with custom subject ID and keep temps

```bash
hippunfold_quick patient_brain.nii.gz /results \
    --subject P001 \
    --keep-temp \
    --cores all
```

### Process with specific hippocampal hemisphere

```bash
hippunfold_quick input.nii.gz /output --hemi R --cores 8
```

## Output

Results are copied to: `{output_dir}/hippunfold/sub-{subject}/`

The structure matches standard hippunfold output:

```
output_dir/
└── hippunfold/
    └── sub-{subject}/
        ├── anat/     # Segmentations and preprocessed images
        ├── surf/     # Surface meshes
        ├── coords/   # Laplace coordinates
        ├── warps/    # Transformation files
        └── qc/       # Quality control images
```

Logs are also copied to: `{output_dir}/logs/`

## Troubleshooting

### "No such file or directory" error
Make sure the input file path is correct and the file exists.

### "Hippunfold failed with exit code 1"
Check the output messages for specific errors. Use `--keep-temp` to inspect the temporary directories.

### Insufficient space in /tmp
If your system's /tmp is small, you can set the `TMPDIR` environment variable:

```bash
export TMPDIR=/path/to/larger/temp
hippunfold_quick input.nii.gz /output
```

## See Also

- [HippUnfold Documentation](https://hippunfold.readthedocs.io/)
- [Installation Guide](installation.md)
- [Docker Guide](docker.md)
