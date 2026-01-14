# Hippunfold Run Guide

This guide describes how to install dependencies and run `hippunfold` using Poetry.

## Installation

### Prerequisites

*   **Python 3.9**: Required due to project constraints (>=3.8,<3.10) and `nibabel` dependency (>=3.9).
*   **Poetry**: Python dependency manager.
*   **Apptainer/Singularity**: Required if using the `--use-singularity` flag (recommended for reproducible environments).

### Steps

1.  **Load Required Modules**:
    ```bash
    module load apptainer
    module load python/3.9
    ```

2.  **Install Poetry** (if not already installed):
    ```bash
    curl -sSL https://install.python-poetry.org | python3.9 -
    ```

3.  **Clone the Repository**:
    ```bash
    git clone https://github.com/khanlab/hippunfold.git
    cd hippunfold
    ```

4.  **Checkout the Correct Branch**:
    ```bash
    git checkout latest-deeplearning
    git branch # Ensure you are on the 'latest-deeplearning' branch
    ```

5.  **Install Dependencies**:
    ```bash
    poetry install
    ```

## Running Hippunfold

To run `hippunfold`, use the following general command structure:

```bash
poetry run hippunfold \
    --modality <modality> \
    <bids_input_dir> \
    <output_dir> \
    participant \
    [options]
```

### Argument Breakdown

*   `--modality <modality>`: Specifies the input modality (e.g., `T1w`, `T2w`).
*   `<bids_input_dir>`: Absolute or relative path to the input BIDS dataset.
*   `<output_dir>`: Path where results will be saved.
*   `participant`: Analysis level (processing participants).
*   `[options]`: Additional flags (see common options below).

### Common Options

*   `--cores <n>`: Number of CPU cores to utilize.
*   `--template <template_name>`: Template to use (e.g., `dHCP` for neonates).
*   `--force_nnunet_model <model_name>`: Force a specific nnUNet model (e.g., `neonateT2w_v2`).
*   `--use-singularity`: Enables execution via Singularity container.
*   `--singularity-args "<args>"`: Arguments passed to Singularity (e.g., `"-B /nfs"` to mount external directories).
    > **Note**: Replace `/nfs` (or whatever path is used) with the most parent directory containing your input data and output locations to avoid Read/Write permission issues or "file not found" errors inside the container.

## Example: Neonate T2w Analysis

Here is a concrete example using test data:

```bash
poetry run hippunfold \
    --modality T2w \
    test_data/bids_dhcp/ \
    test_output/bids_dhcp_nnunet_T2w/ \
    participant \
    --cores 8 \
    --template dHCP \
    --force_nnunet_model neonateT2w_v2 \
    --use-singularity \
    --singularity-args "-B /nfs"
```
