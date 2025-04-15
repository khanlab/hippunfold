#!/usr/bin/env python3
import shutil
import subprocess
import sys
import tempfile
import os
from argparse import ArgumentParser
from pathlib import Path


def check_conda_installation():
    try:
        conda_version = subprocess.run(
            ["conda", "--version"], capture_output=True, text=True, check=True
        )
        print(f"Conda is installed (version: {conda_version.stdout.strip()})")
        return True
    except FileNotFoundError:
        error_message = "Conda is not installed on your system or not found in PATH."
        print(error_message)
        return False


def gen_parser():
    parser = ArgumentParser(description="Run hippunfold for a single subject.")

    # command line arguments for hippunfold-quick
    parser.add_argument(
        "-i","--input", required=True, help="File path to your input NIfTI image. Must have .nii.gz extension."
    )
    parser.add_argument(
        "-o","--output", required=True, help="Path to your desired output folder."
    )
    parser.add_argument("-s","--subject", required=True, help="Subject ID (e.g., 001).")
    parser.add_argument(
        "-m", "--modality",
        required=True,
        choices=["T1w", "T2w"],  #currently hardcoded to put things in anat
        help="Image modality.",
    )
    parser.add_argument(
        "--temp-dir",
        default=None,
        help="Optional temporary directory. If not specified, a system temp directory will be used.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Execute a dry run without actually running the full pipeline.",
    )

    return parser


def main():

    script_path = Path(__file__).resolve().parent / "run.py"

    if not check_conda_installation():
        print("Please activate conda and install snakebids to continue using hippunfold-quick.")
        sys.exit(1)

    if "SNAKEMAKE_PROFILE" in os.environ:
        del os.environ["SNAKEMAKE_PROFILE"]
    
    args = gen_parser().parse_args()

    # set temp dir if specified, else use python tempfile
    if args.temp_dir:
        prefix = args.temp_dir
    else:
        prefix = None

    with tempfile.TemporaryDirectory(prefix=prefix) as temp_dir:

        # create subject folder
        subject_folder = Path(temp_dir) / "anat" / f"sub-{args.subject}"
        subject_folder.mkdir(parents=True,exist_ok=True)

        # create new file name
        input_filename = f"sub-{args.subject}_{args.modality}.nii.gz"

        # add file name to the created subject folder
        temp_input_file = subject_folder / input_filename

        # copy the input file
        shutil.copy(args.input, temp_input_file)

        # run hippunfold
        command = [
            script_path,
            temp_dir,
            args.output,
            "participant",
            "-c",
            "all",
            "--force-output",
            "--nolock",
            "--modality",
            args.modality,
            "--use-conda",
            #"--quiet",
        ]

        if args.dry_run:
            command.append("-n")

        # run the command
        try:
            subprocess.run(command, check=True)
            print("hippunfold completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    main()
