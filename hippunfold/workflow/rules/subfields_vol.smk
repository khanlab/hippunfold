
#get subfield labels using volumetric coords:
rule label_subfields_from_vol_coords_corobl:
    input:  
        subfields_mat = os.path.join(config['snakemake_dir'],'hippocampal_autotop','misc','BigBrain_ManualSubfieldsUnfolded.mat'),
        nii_ap = bids(root='work',**config['subj_wildcards'],suffix='autotop/coords-AP.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        nii_pd = bids(root='work',**config['subj_wildcards'],suffix='autotop/coords-PD.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    params:
        mat_name = 'subfields_avg' #avg bigbrain over L/R hemis
    output:
        nii_label = bids(root='work',datatype='seg_{modality}',desc='subfields',suffix='dseg.nii.gz', space='corobl',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    script: '../scripts/label_subfields_from_vol_coords.py'


#-- resampling to T1w native space
rule resample_subfields_to_T1w:
    input:
        nii = bids(root='work',datatype='seg_{modality}',desc='subfields',suffix='dseg.nii.gz', space='corobl',from_='volume',hemi='{hemi}', **config['subj_wildcards']),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz')
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 
 
rule resample_postproc_to_T1w:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz')
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='postproc',space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 
  
rule resample_unet_to_T1w:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/manual_lbl.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='T1w.nii.gz')
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='unet',space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 


#create ref space for hires crop in native space
# TODO:  expose the resampling factor and size as cmd line args
rule create_native_crop_ref:
    input:
        seg = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    params:
        resample = '400%',
        pad_to = '192x256x256vox'
    output:
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'c3d {input} -binarize -interpolation NearestNeighbor -trim 0vox -resample {params.resample} -pad-to {params.pad_to} 0 {output}'
  
 

rule resample_unet_native_crop:
    """ This is either nnUnet or manual seg, since currently nnUnet is passed to autotop with the manual seg param"""
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/manual_lbl.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='unet',space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 

     
rule resample_postproc_native_crop:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='postproc',space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 


rule resample_subfields_native_crop:
    input:
        nii = bids(root='work',datatype='seg_{modality}',desc='subfields',suffix='dseg.nii.gz', space='corobl',from_='volume',hemi='{hemi}', **config['subj_wildcards']),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='cropT1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation MultiLabel -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 


rule resample_coords_native_crop:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/coords-{dir}.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz', space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.nii} -o {output.nii} -r {input.ref}  -t [{input.xfm},1]' 


#add srlm, cyst, dg from postproc labels to subfields
#input dg label 8, output 6
#input srlm label 2, output 7
#input cyst label 7, output 8

#first remap tissue labels to get three sep labels
# then, we just need to add those in, using max(old,new) to override old with new in case of conflict
rule combine_tissue_subfield_labels_native_crop:
    input:
        tissue = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='postproc',space='cropT1w',hemi='{hemi}', **config['subj_wildcards']),
        subfields = bids(root='work',datatype='seg_{modality}',desc='subfields',suffix='dseg.nii.gz', space='cropT1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    params:
        remap_dg = '-threshold 8 8 6 0 -popas dg',
        remap_srlm = '-threshold 2 2 7 0 -popas srlm',
        remap_cyst = '-threshold 7 7 8 0 -popas cyst',
    output:
        combined = bids(root='work',datatype='seg_{modality}',desc='subfieldswithtissue',suffix='dseg.nii.gz', space='cropT1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 
        'c3d {input.tissue} -dup {params.remap_dg} -dup {params.remap_srlm} {params.remap_cyst} {input.subfields} -push dg -max -push srlm -max -push cyst -max -o {output}'
        
rule combine_tissue_subfield_labels_native:
    input:
        tissue = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='postproc',space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        subfields = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    params:
        remap_dg = '-threshold 8 8 6 0 -popas dg',
        remap_srlm = '-threshold 2 2 7 0 -popas srlm',
        remap_cyst = '-threshold 7 7 8 0 -popas cyst',
    output:
        combined = bids(root='work',datatype='seg_{modality}',desc='subfieldswithtissue',suffix='dseg.nii.gz', space='T1w',from_='volume',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 
        'c3d {input.tissue} -dup {params.remap_dg} -dup {params.remap_srlm} {params.remap_cyst} {input.subfields} -push dg -max -push srlm -max -push cyst -max -o {output}'
        

"""
#right now this uses same labels for each, need to change this to a new lut
rule combine_lr_subfields:
    input:
        left = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='L', **config['subj_wildcards']),
        right = bids(root='work',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='R', **config['subj_wildcards'])
    output:
        combined = bids(root='results',datatype='seg_{modality}',suffix='dseg.nii.gz', desc='subfields',space='T1w', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 'c3d {input} -add -o {output}'
 """




