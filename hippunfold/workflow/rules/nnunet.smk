
rule extract_model:
    input: 
        model_tar = lambda wildcards: config['nnunet_model'][wildcards.modality]
    output: 
        model_dir = directory('nnunet_models/{modality}'),
    shell: 'mkdir -p {output} && tar -xvf {input} -C {output}'

def get_nnunet_input (wildcards):
    if wildcards.modality == 'T2w':
        nii = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T2w.nii.gz',desc='cropped',space='corobl',hemi='{hemi}'),
    else:
        raise ValueError('modality not supported for nnunet!')
    return nii


rule copy_img_to_temp:
    input: get_nnunet_input
    params:
        temp_img_dir = bids(root='work',suffix='tempimg',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}')
    output: 
        temp_img = bids(root='work',suffix='tempimg/temp_0000.nii.gz',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}')
    threads: 32
    group: 'subj'
    shell: 'mkdir -p {params.temp_img_dir} &&  cp {input} {output.temp_img}'

def parse_task_from_tar (wildcards, input):
    import re
    match = re.search('Task[0-9]{3}_[a-zA-Z0-9]+',input.model_tar)
    if match:
        task = match.group(0)
    else:
        raise ValueError('cannot parse Task from model tar')
    return task



def parse_chkpnt_from_tar (wildcards, input):
    import re
    match = re.search('^.*\.(\w+)\.tar',input.model_tar)
    if match:
        chkpnt = match.group(1)
    else:
        raise ValueError('cannot parse chkpnt from model tar')
    return chkpnt

rule run_inference:
    input: 
        in_img = bids(root='work',suffix='tempimg/temp_0000.nii.gz',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}'),
        model_dir = 'nnunet_models/{modality}',
        model_tar = lambda wildcards: config['nnunet_model'][wildcards.modality]
    params:
        in_folder = bids(root='work',suffix='tempimg',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}'),
        out_folder = bids(root='work',suffix='templbl',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}'),
        task = parse_task_from_tar,
        chkpnt = parse_chkpnt_from_tar,
    output: 
        tmp_lbl = bids(root='work',suffix='templbl/temp.nii.gz',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}')
    threads: 8 
    resources:
        gpus = 1,
        mem_mb = 64000,
        time = 30,
    group: 'subj'
    shell: 'export RESULTS_FOLDER={input.model_dir} &&'
           'export nnUNet_n_proc_DA={threads} &&'
           'nnUNet_predict -i {params.in_folder} -o {params.out_folder} -t {params.task} -chk {params.chkpnt}' # --disable_tta'

rule copy_temp_to_lbl:
    input:
        tmp_lbl = bids(root='work',suffix='templbl/temp.nii.gz',hemi='{hemi}',**config['subj_wildcards'],modality='{modality}')
    output:
        nnunet_seg = bids(root='work',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='nnunet',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')
    threads: 32
    group: 'subj'
    shell: 'cp {input} {output}'


