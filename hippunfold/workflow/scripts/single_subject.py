from argparse import ArgumentParser
import os
import tempfile
import shutil
import subprocess

def gen_parser():
    parser = ArgumentParser(description="Run hippunfold for a single subject.")

    # command line arguments for hippunfold-singlesubject
    parser.add_argument(
        "--input", 
        required=True, 
        help="File path to your input NIfTI image."
    )
    parser.add_argument(
        "--output", 
        required=True, 
        help="Path to your desired output folder."
    )
    parser.add_argument(
        "--subject", 
        required=True, 
        help="Subject ID (e.g., 001)."
    )
    parser.add_argument(
        "--modality", 
        required=True, 
        choices=["T1w", "T2w", "hippb500", "dsegtissue"], 
        help="Image modality."
    )
    parser.add_argument(
        "--temp-dir", 
        default=None, 
        help="Optional temporary directory. If not specified, a system temp directory will be used."
    )
    parser.add_argument(
        "--use-singularity",
        action="store_true",
        help="If specified, will run hippunfold using Singularity. "
            "This requires Singularity to be installed on your system."
    )
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help="If specified, will run hippunfold using Conda environments. "
            "This requires Conda to be installed on your system."
    )
    parser.add_argument(
        "-np", "--dry-run",
        action="store_true",
        help="Execute a dry run without actually running the full pipeline."
    )

    return parser

def main():
    args = gen_parser().parse_args()
     
    # set temp dir if specified, else use python tempfile
    if args.temp_dir:
        temp_dir = args.temp_dir
        os.makedirs(temp_dir, exist_ok=True)
    else:
        temp_dir = tempfile.mkdtemp()
    
    # create subject folder 
    subject_folder = os.path.join(temp_dir, f"anat/sub-{args.subject}")
    os.makedirs(subject_folder, exist_ok=True) 

    # create new file name
    input_filename = f"sub-{args.subject}_{args.modality}.nii.gz"
    
    # add file name to the created subject folder
    temp_input_file = os.path.join(subject_folder, input_filename)
    
    # copy the input file
    shutil.copy(args.input, temp_input_file) 
     
    # run hippunfold
    command = [
        "hippunfold",
        temp_dir,
        args.output,
        "participant",
        "-c", "all",
        "--force-output",
        "--nolock",
        "--modality", args.modality
    ]

    if args.use_singularity:
        command.append("--use-singularity")
    if args.use_conda:
        command.append("--use-conda")
    if args.dry_run:
        command.append("-np")

    # run the command 
    try:
        subprocess.run(command, check=True)
        print("hippunfold completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    finally:
        # delete the temp dir if not user specified
        if not args.temp_dir:
            shutil.rmtree(temp_dir)
 
if __name__ == "__main__":
    main()