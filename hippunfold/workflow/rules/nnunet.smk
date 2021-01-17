import re

def get_nnunet_input (wildcards):
    if wildcards.modality == 'T2w':
        nii = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T2w.nii.gz',desc='cropped',space='corobl',hemi='{hemi}'),
    else:
        raise ValueError('modality not supported for nnunet!')
    return nii


def parse_task_from_tar (wildcards, input):
    match = re.search('Task[0-9]{3}_[a-zA-Z0-9]+',input.model_tar)
    if match:
        task = match.group(0)
    else:
        raise ValueError('cannot parse Task from model tar')
    return task

def parse_chkpnt_from_tar (wildcards, input):
    match = re.search('^.*\.(\w+)\.tar',input.model_tar)
    if match:
        chkpnt = match.group(1)
    else:
        raise ValueError('cannot parse chkpnt from model tar')
    return chkpnt

rule run_inference:
    """ This rule REQUIRES a GPU -- will need to modify nnUnet code to create an alternate for CPU-based inference
        It also runs in an isolated folder (shadow), with symlinks to inputs in that folder, copying over outputs once complete, so temp files are not retained"""
    input: 
        in_img = get_nnunet_input,
        model_tar = lambda wildcards: config['nnunet_model'][wildcards.modality]
    params:
        temp_img = 'tempimg/temp_0000.nii.gz',
        temp_lbl = 'templbl/temp.nii.gz',
        model_dir = 'tempmodel',
        in_folder = 'tempimg',
        out_folder = 'templbl',
        task = parse_task_from_tar,
        chkpnt = parse_chkpnt_from_tar,
    output: 
        nnunet_seg = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='nnunet',space='corobl',hemi='{hemi,Lflip|R}')
    shadow: 'minimal' 
    threads: 16 
    resources:
        gpus = 1,
        mem_mb = 32000,
        time = 30,
    group: 'subj'
    shell: 'mkdir -p {params.model_dir} {params.in_folder} {params.out_folder} && ' #create temp folders
           'cp -v {input.in_img} {params.temp_img} && ' #cp input image to temp folder
           'tar -xvf {input.model_tar} -C {params.model_dir} && ' #extract model
           'export RESULTS_FOLDER={params.model_dir} && ' #set nnunet env var to point to model
           'export nnUNet_n_proc_DA={threads} && ' #set threads
           'nnUNet_predict -i {params.in_folder} -o {params.out_folder} -t {params.task} -chk {params.chkpnt} && ' # run inference
           'cp -v {params.temp_lbl} {output.nnunet_seg}' #copy from temp output folder to final output

