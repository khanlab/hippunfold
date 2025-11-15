# HippUnfold Development Setup Guide

This guide walks you through setting up HippUnfold for neonatal brain segmentation using Pixi.

## Prerequisites

- Git
- Linux/macOS system
- ~10GB disk space for dependencies
- **Internet connection** (required for model downloads and atlas sync)

## Step 1: Clone the Repository

```bash
git clone https://github.com/khanlab/hippunfold.git
cd hippunfold
```

## Step 2: Checkout Development Branch

```bash
git checkout dev-deeplearning
```

Verify you're on the correct branch:
```bash
git branch
# Should show:
# * dev-deeplearning
#   master
```

## Step 3: Install Pixi

Pixi is a fast, cross-platform package manager that will handle all dependencies.

### Download and Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### Activate Pixi in Current Shell

```bash
export PATH="$HOME/.pixi/bin:$PATH"
```

To make this permanent, add to your `~/.bashrc` or `~/.bash_profile`:
```bash
echo 'export PATH="$HOME/.pixi/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Verify Installation

```bash
pixi --version
```

## Step 4: Set Up HippUnfold Environment

Install all required dependencies:

```bash
pixi install
```

This will download and install all required dependencies in an isolated environment.

## Step 5: Running HippUnfold for Neonatal Data

You can run HippUnfold in two ways:

### Option 1: Using `pixi run` (Recommended)

Run commands directly without entering a shell:

```bash
pixi run hippunfold --modality <modality> <input_bids_dir> <output_dir> participant --cores <N> --template dHCP --force_nnunet_model <model>
```

### Option 2: Using `pixi shell`

Enter the development environment shell, then run commands:

```bash
pixi shell
hippunfold --modality <modality> <input_bids_dir> <output_dir> participant --cores <N> --template dHCP --force_nnunet_model <model>
exit  # when done
```

## Step 6: Neonatal Models

### Neonatal T1w Model (nnUNet v2)

```bash
pixi run hippunfold --modality T1w <input_dir> <output_dir> participant \
  --cores 10 \
  --force_nnunet_model neonateT1w_v2 \
  --template dHCP
```

### Neonatal T2w Model (nnUNet v2)

```bash
pixi run hippunfold --modality T2w <input_dir> <output_dir> participant \
  --cores 10 \
  --force_nnunet_model neonateT2w_v2 \
  --template dHCP
```

## Step 7: Example with Test Data

### Neonatal T2w with test data

```bash
pixi run hippunfold --modality T2w test_data/bids_dhcp/ output/dhcp_T2w participant \
  --cores 10 \
  --force_nnunet_model neonateT2w_v2 \
  --template dHCP
```

### Neonatal T1w with test data

```bash
pixi run hippunfold --modality T1w test_data/bids_dhcp/ output/dhcp_T1w participant \
  --cores 10 \
  --force_nnunet_model neonateT1w_v2 \
  --template dHCP
```

### GPU Support

To enable GPU acceleration (if available):

```bash
pixi run hippunfold --modality T2w <input_dir> <output_dir> participant \
  --cores 10 \
  --use_gpu \
  --force_nnunet_model neonateT2w_v2 \
  --template dHCP
```

## Daily Workflow

After initial setup, your daily workflow is:

### Quick Run (pixi run)
```bash
cd /path/to/hippunfold
pixi run hippunfold --modality <modality> <input> <output> participant \
  --cores <N> \
  --force_nnunet_model <model> \
  --template dHCP
```

### Interactive Shell (pixi shell)
```bash
cd /path/to/hippunfold
pixi shell
hippunfold --modality <modality> <input> <output> participant \
  --cores <N> \
  --force_nnunet_model <model> \
  --template dHCP
exit  # when done
```
