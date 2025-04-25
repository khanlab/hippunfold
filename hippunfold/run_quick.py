#!/usr/bin/env python3
import os
import shutil
import subprocess
import sys
import tempfile
from argparse import ArgumentParser
from pathlib import Path

IMAGE_MODALITY = {
    "T1w": {
        "suffix": "T1w",
        "use_derivatives": False,
    },
    "T2w": {
        "suffix": "T2w",
        "use_derivatives": False,
    },
    "dsegtissue": {
        "suffix": "dseg",
        "use_derivatives": True,
    },
}


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
        "-i",
        "--input",
        required=True,
        help="File path to your input NIfTI image. Must have .nii.gz extension.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to your desired output folder."
    )
    parser.add_argument(
        "-s", "--subject", required=True, help="Subject ID (e.g., 001)."
    )
    parser.add_argument(
        "-m",
        "--modality",
        required=True,
        choices=list(
            IMAGE_MODALITY.keys()
        ),  # currently hardcoded to put things in anat
        help="Image modality - choose between: " + ", ".join(IMAGE_MODALITY.keys()),
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
    parser.add_argument(
        "--hemi",
        required=False,
        choices=["L", "R"],
        help="Brain hemisphere, tequired for dsegtissue modality. Choose between 'L' or 'R'.",
    )

    return parser


def main():

    script_path = Path(__file__).resolve().parent / "run.py"

    if not check_conda_installation():
        print(
            "Please activate conda and install snakebids to continue using hippunfold-quick."
        )
        sys.exit(1)

    if "SNAKEMAKE_PROFILE" in os.environ:
        del os.environ["SNAKEMAKE_PROFILE"]

    args = gen_parser().parse_args()

    # if the user selects dsegtissue but doesn't specify hemi, throw an error
    if args.modality == "dsegtissue" and not args.hemi:
        print(
            "Error: The 'hemi' argument is required when using the 'dsegtissue' modality."
        )
        sys.exit(1)

    # set temp dir if specified, else use python tempfile
    if args.temp_dir:
        prefix = args.temp_dir
    else:
        prefix = None

    with tempfile.TemporaryDirectory(prefix=prefix) as temp_dir:

        # create subject folder
        subject_folder = Path(temp_dir) / "anat" / f"sub-{args.subject}"
        subject_folder.mkdir(parents=True, exist_ok=True)

        # create new file name
        # if dsegtissue, desc-tissue needs to be added for bids format
        if args.modality == "dsegtissue":
            input_filename = f"sub-{args.subject}_hemi-{args.hemi}_desc-tissue_{IMAGE_MODALITY[args.modality]['suffix']}.nii.gz"
        else:
            input_filename = (
                f"sub-{args.subject}_{IMAGE_MODALITY[args.modality]['suffix']}.nii.gz"
            )

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
            "--quiet",
        ]

        if args.dry_run:
            command.append("-n")

        if IMAGE_MODALITY[args.modality]["use_derivatives"]:
            # need to have desc file in bids dir to use --derivatives
            desc_file = Path(temp_dir) / "dataset_description.json"
            desc_file.write_text(
                '{"Name": "Generated Derivatives", '
                '"BIDSVersion": "1.0.2", '
                '"GeneratedBy": [{"Name": "hippunfold-quick"}]}'
            )
            command.append("--derivatives")
            command.append(Path(temp_dir))
            command.append("hemi")
            command.append(args.hemi)

        # run the command
        try:
            subprocess.run(command, check=True)
            print("hippunfold completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")


if __name__ == "__main__":
    main()
