#!/usr/bin/env python3
"""Quick run script for hippunfold that uses tempdir for improved efficiency.

This script creates a temporary BIDS dataset from input files and runs hippunfold
using temporary directories for both input and output. After successful completion,
results are copied back to the specified output directory. This is particularly
useful when working with network drives where local temporary storage provides
better I/O performance.

Usage:
    hippunfold_quick INPUT_FILE OUTPUT_DIR [OPTIONS]

    Basic example:
        hippunfold_quick sub-001_T2w.nii.gz /path/to/output

    With custom subject ID and modality:
        hippunfold_quick brain.nii.gz /output --subject 123 --modality T1w

    With additional hippunfold options:
        hippunfold_quick input.nii.gz /output -n --cores 8 --hemi L

    Keep temporary files for debugging:
        hippunfold_quick input.nii.gz /output --keep-temp

The script performs the following steps:
    1. Creates a temporary BIDS dataset from the input file
    2. Creates a temporary output directory in the same temp location
    3. Runs hippunfold using the temporary directories
    4. Copies the subject results from temp to the final output directory
    5. Cleans up temporary files (unless --keep-temp is specified)

This approach improves efficiency when:
    - Output directory is on a network drive (e.g., NFS, CIFS)
    - You want to isolate the workflow from the final output location
    - You need to process single images without manually creating BIDS structure
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def create_bids_dataset(input_file, subject_id, tempdir, modality="T2w"):
    """Create a minimal BIDS dataset in tempdir from a single input file.

    Args:
        input_file: Path to the input image file
        subject_id: Subject identifier
        tempdir: Path to temporary directory for BIDS dataset
        modality: Image modality (T1w, T2w, etc.)

    Returns:
        Path to the created BIDS directory
    """
    bids_dir = Path(tempdir) / "bids_input"
    bids_dir.mkdir(parents=True, exist_ok=True)

    # Create dataset_description.json
    dataset_desc = bids_dir / "dataset_description.json"
    dataset_desc.write_text(
        '{"Name": "Quick hippunfold run", '
        '"BIDSVersion": "1.6.0", '
        '"DatasetType": "raw"}\n'
    )

    # Create subject directory
    subj_dir = bids_dir / f"sub-{subject_id}" / "anat"
    subj_dir.mkdir(parents=True, exist_ok=True)

    # Copy or link input file
    output_filename = f"sub-{subject_id}_{modality}.nii.gz"
    output_path = subj_dir / output_filename
    shutil.copy2(input_file, output_path)

    # Create minimal JSON sidecar for BIDS compliance
    json_filename = f"sub-{subject_id}_{modality}.json"
    json_path = subj_dir / json_filename
    json_path.write_text('{"Modality": "' + modality + '"}\n')

    return bids_dir


def main():
    """Main function for quick hippunfold runs."""
    parser = argparse.ArgumentParser(
        description="Quick hippunfold run using temporary directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        help="Path to input image file (NIfTI format)",
    )
    parser.add_argument(
        "output_dir",
        help="Path to final output directory",
    )
    parser.add_argument(
        "--subject",
        default="001",
        help="Subject identifier (default: 001)",
    )
    parser.add_argument(
        "--modality",
        choices=["T1w", "T2w", "hippb500"],
        default="T2w",
        help="Image modality (default: T2w)",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary directories after completion (for debugging)",
    )

    # Parse known args and pass rest to hippunfold
    args, hippunfold_args = parser.parse_known_args()

    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file does not exist: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    final_output_dir = Path(args.output_dir).absolute()
    final_output_dir.mkdir(parents=True, exist_ok=True)

    # Create temporary directory for workflow
    temp_base_dir = tempfile.mkdtemp(prefix="hippunfold_quick_")

    try:
        print(f"Using temporary directory: {temp_base_dir}")

        # Create BIDS dataset in tempdir
        print(f"Creating temporary BIDS dataset for subject {args.subject}...")
        bids_dir = create_bids_dataset(
            args.input_file, args.subject, temp_base_dir, args.modality
        )

        # Create output directory in tempdir
        temp_output_dir = Path(temp_base_dir) / "hippunfold_output"
        temp_output_dir.mkdir(parents=True, exist_ok=True)

        print("Running hippunfold...")
        print(f"  Input (temp):  {bids_dir}")
        print(f"  Output (temp): {temp_output_dir}")

        # Prepare arguments for hippunfold
        hippunfold_cli_args = [
            str(bids_dir),
            str(temp_output_dir),
            "participant",
            "--modality",
            args.modality,
            "--participant_label",
            args.subject,
        ] + hippunfold_args

        # Run hippunfold using subprocess
        # Build the command to run hippunfold
        # Use the same Python interpreter and call run.py directly
        run_py_path = os.path.join(os.path.dirname(__file__), "run.py")
        cmd = [sys.executable, run_py_path] + hippunfold_cli_args

        result = subprocess.run(cmd, capture_output=False)

        if result.returncode != 0:
            raise RuntimeError(f"Hippunfold failed with exit code {result.returncode}")

        # Copy results from temp output to final output
        print(f"\nCopying results to final output directory: {final_output_dir}")

        # Copy the hippunfold results directory
        temp_hippunfold_dir = temp_output_dir / "hippunfold"
        if temp_hippunfold_dir.exists():
            final_hippunfold_dir = final_output_dir / "hippunfold"
            final_hippunfold_dir.mkdir(parents=True, exist_ok=True)

            # Copy only the subject directory
            subj_dir_name = f"sub-{args.subject}"
            temp_subj_dir = temp_hippunfold_dir / subj_dir_name

            if temp_subj_dir.exists():
                final_subj_dir = final_hippunfold_dir / subj_dir_name
                if final_subj_dir.exists():
                    print(f"  Removing existing subject directory: {final_subj_dir}")
                    shutil.rmtree(final_subj_dir)

                print(f"  Copying {subj_dir_name}...")
                shutil.copytree(temp_subj_dir, final_subj_dir)
            else:
                print(
                    "Warning: Subject directory not found in temp output: "
                    f"{temp_subj_dir}",
                    file=sys.stderr,
                )
        else:
            print(
                "Warning: hippunfold directory not found in temp output: "
                f"{temp_hippunfold_dir}",
                file=sys.stderr,
            )

        # Optionally copy logs
        temp_logs_dir = temp_output_dir / "logs"
        if temp_logs_dir.exists():
            final_logs_dir = final_output_dir / "logs"
            final_logs_dir.mkdir(parents=True, exist_ok=True)

            # Copy subject-specific logs
            for log_file in temp_logs_dir.glob(f"*sub-{args.subject}*"):
                shutil.copy2(log_file, final_logs_dir)

        print("\nHippunfold quick run completed successfully!")
        print(f"Results are in: {final_output_dir}")

    except Exception as e:
        print(f"\nError during hippunfold execution: {e}", file=sys.stderr)
        if args.keep_temp:
            print(
                f"Temporary directory preserved at: {temp_base_dir}",
                file=sys.stderr,
            )
        raise
    finally:
        # Clean up temp directory unless --keep-temp is specified
        if not args.keep_temp:
            print("\nCleaning up temporary directory...")
            shutil.rmtree(temp_base_dir, ignore_errors=True)
        else:
            print(f"\nTemporary directory preserved at: {temp_base_dir}")


if __name__ == "__main__":
    main()
