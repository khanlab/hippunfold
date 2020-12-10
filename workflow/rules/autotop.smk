import os
    

#creates command based on mcr or matlab
def get_autotop_cmd (wildcards, input, output):
    autotop_dir = os.path.join(config['snakemake_dir'],'hippocampal_autotop')
    singularity_cmd = f"singularity exec -B {autotop_dir}:/src -e {config['singularity']['autotop']}" 

    if config['use_mcr'] == True:
        cmd = f"{singularity_cmd} /src/mcr_v97/run_AutoTops_TransformAndRollOut.sh /opt/mcr/v97 "\
                f"{input.nii} {output.out_dir} '' {config['cnn_model'][wildcards.modality]}"
    else:
        set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
        set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"
        set_autotop_dir = f"SINGULARITYENV_AUTOTOP_DIR={autotop_dir}"

        cmd = f"{set_matlab_lic} {set_java_home} {set_autotop_dir} {singularity_cmd} "\
                f"{config['matlab_bin']} -batch \"addpath(genpath('{autotop_dir}')); "\
                f"AutoTops_TransformAndRollOut('{input.nii}','{output.out_dir}','{config['cnn_model'][wildcards.modality]}')\""
    return cmd   



def get_autotop_input (wildcards):
    if wildcards.modality == 'T2w':
        nii = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T2w.nii.gz',desc='cropped',space='corobl',hemi='{hemi}'),
    elif wildcards.modality == 'T1w':
        nii = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='InvT1w.nii.gz',desc='cropped',space='corobl',hemi='{hemi}'),
    elif wildcards.modality == 'b500':
        nii = bids(root='work/preproc_dwi',suffix='b500.nii.gz',desc='cropped',datatype='dwi',**config['subj_wildcards'],space='corobl',hemi='{hemi}'),
    elif wildcards.modality == 'b1000':
        nii = bids(root='work/preproc_dwi',suffix='b1000.nii.gz',desc='cropped',datatype='dwi',**config['subj_wildcards'],space='corobl',hemi='{hemi}'),
    else:
        nii = ''

    return nii



rule run_autotop:
    input:
        nii = get_autotop_input
    params:
        autotop_cmd = get_autotop_cmd
    output:
        out_dir = directory(bids(root='work',**config['subj_wildcards'],suffix='autotop',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')),
        subfields = bids(root='work',**config['subj_wildcards'],suffix='autotop/subfields-BigBrain.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2native_extrap = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native_extrapolateNearest.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        gii = expand(bids(root='work',suffix='autotop/{surfname}.unfoldedtemplate.surf.gii',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),surfname=['inner','outer','midthickness'],allow_missing=True),
        coords = expand(bids(root='work',suffix='autotop/coords-{dir}.nii.gz',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),dir=['AP','PD','IO'],allow_missing=True)
    threads: 8
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='autotop.txt')
    shell:
        'SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        "{params.autotop_cmd} &> {log}"

def get_autotop_inputseg_cmd (wildcards, input, output):
    autotop_dir = os.path.join(config['snakemake_dir'],'hippocampal_autotop')
    singularity_cmd = f"singularity exec -B {autotop_dir}:/src -e {config['singularity']['autotop']}" 

    if config['use_mcr'] == True:
        cmd = f"{singularity_cmd} /src/mcr_v97/run_AutoTops_TransformAndRollOut.sh /opt/mcr/v97 "\
                f"{input.nii} {output.out_dir} '{input.seg}' {config['cnn_model'][wildcards.modality]}"
    else:
        set_matlab_lic = f"SINGULARITYENV_MLM_LICENSE_FILE={config['mlm_license_file']}"
        set_java_home = f"SINGULARITYENV_JAVA_HOME={config['java_home']}"
        set_autotop_dir = f"SINGULARITYENV_AUTOTOP_DIR={autotop_dir}"

        cmd = f"{set_matlab_lic} {set_java_home} {set_autotop_dir} {singularity_cmd} "\
                f"{config['matlab_bin']} -batch \"addpath(genpath('{autotop_dir}')); "\
                f"AutoTops_TransformAndRollOut('{input.nii}','{output.out_dir}','{config['cnn_model'][wildcards.modality]},'{input.seg}'')\""
    return cmd   

rule run_autotop_inputseg:
    input:
        nii = get_autotop_input,
        seg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',from_='{modality}'),
    params:
        autotop_cmd = get_autotop_inputseg_cmd
    output:
        out_dir = directory(bids(root='work',**config['subj_wildcards'],suffix='autotop',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}')),
        subfields = bids(root='work',**config['subj_wildcards'],suffix='autotop/subfields-BigBrain.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}'),
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}'),
        gii = expand(bids(root='work',suffix='autotop/{surfname}.unfoldedtemplate.surf.gii',desc='cropped', space='corobl',hemi='{{hemi}}',modality='seg{{modality}}', **config['subj_wildcards']),surfname=['inner','outer','midthickness'],allow_missing=True)
    threads: 8
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}',suffix='autotop.txt')
    shell:
        'SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        "{params.autotop_cmd} &> {log}"



#full-grid correction of unfolded space
rule map_to_full_grid:
    input: 
        autotop_dir = bids(root='work',**config['subj_wildcards'],suffix='autotop',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    params:
        script = os.path.join(config['snakemake_dir'],'hippocampal_autotop','tools','warps_gifti','mapUnfoldToFullGrid.sh')
    output:
        warp_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        '{params.script} {input.autotop_dir} {input.autotop_dir}'


#rule to unflip a nifti
rule unflip_autotop_nii:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/{filename}.nii.gz',desc='cropped',space='corobl',hemi='{hemi}flip',modality='{modality}')
    output:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/{filename}.nii.gz',desc='cropped',space='corobl',hemi='{hemi,L}',modality='{modality}')
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell: 'c3d {input} -flip x {output}'


   

rule resample_subfields_to_T1w:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/subfields-BigBrain.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz')
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 

#create ref space for hires crop in native space
# TODO:  expose the resampling factor and size as cmd line args
rule create_native_crop_ref:
    input:
        seg = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'c3d {input} -binarize -interpolation NearestNeighbor -trim 0vox -resample 400% -pad-to 160x320x160vox 0 {output}'
  
rule resample_subfields_native_crop:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/subfields-BigBrain.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 


rule resample_coords_native_crop:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/coords-{dir}.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz', space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 




#right now this uses same labels for each, need to change this to a new lut
rule combine_lr_subfields:
    input:
        left = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='L', **config['subj_wildcards']),
        right = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='R', **config['subj_wildcards'])
    output:
        combined = bids(root='results',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w', **config['subj_wildcards'])
    container: config['singularity']['prepdwi']
    group: 'subj'
    shell: 'c3d {input} -add -o {output}'
 



