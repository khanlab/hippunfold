rule create_native_crop_ref:
    """Create ref space for hires crop in native space
        TODO:  expose the resampling factor and size as cmd line args"""
    input:
        seg = bids(root='results',datatype='seg',suffix='dseg.nii.gz', desc='subfields',space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    params:
        resample = '400%',
        pad_to = '192x256x256vox'
    output:
        ref = bids(root='work',datatype='seg',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'c3d {input} -binarize -interpolation NearestNeighbor -trim 0vox -resample {params.resample} -pad-to {params.pad_to} 0 -type uchar {output}'
  

rule compose_warps_corobl2unfold_rhemi:
    """ Compose corobl to unfold (unfold-template), for right hemi (ie no flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
    output:
        corobl2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.native2unfold}'


rule compose_warps_corobl2unfold_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt')
    output:
        corobl2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.native2unfold} {input.flipLR_xfm}'




rule compose_warps_unfold2corobl_rhemi:
    """ Compose unfold (unfold-template) to corobl, for right hemi (ie no flip)"""
    input:
        unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native_RPI.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        unfold_ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        ref = bids(root='work',datatype='seg',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='R', **config['subj_wildcards']),
    output:
        unfold2corobl = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='unfold',to='corobl',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2corobl},1] -r {input.unfold2native} -t {input.unfold2native} -i {input.unfold_ref} -v'


rule compose_warps_unfold2corobl_lhemi:
    """ Compose unfold (unfold-template) to corobl, for left hemi (ie with flip)"""
    input:
        unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native_RPI.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt'),
        unfold_ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        ref = bids(root='work',datatype='seg',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='L', **config['subj_wildcards']),
    output:
        unfold2corobl = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='unfold',to='corobl',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2corobl},1] -r {input.ref} -t {input.flipLR_xfm} -t {input.unfold2native} -i {input.unfold_ref} -v'


 
rule compose_warps_t1_to_unfold:
    """ Compose warps from T1w to unfold """
    input:
        corobl2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='corobl',to='unfold',mode='image'),
        ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
    output: 
        bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='T1w',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output} -R {input.ref} {input.corobl2unfold} {input.t1w2corobl}'

print(bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='cropT1w',mode='image'))

rule compose_warps_unfold_to_cropt1:
    """ Compose warps from unfold to cropT1w """
    input:
        unfold2corobl = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='corobl',mode='image'),
        ref = bids(root='work',datatype='seg',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        unfold_ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
    output: 
        unfold2cropt1 = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='T1w',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2cropt1},1] -r {input.ref} -t [{input.t1w2corobl},1] -t {input.unfold2corobl} -i {input.unfold_ref} -v'




rule fix_unfold2native_orient:
    """Fix orientation issue with WarpITK_unfold2native.nii.gz (LPS to RPI)
        TODO: fix this at the source (autotop)"""
    input:
        bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    output:
        bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native_RPI.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    group: 'subj'
    shell: 'c3d -mcs {input} -orient RPI -omc {output}'


rule resample_T1w_to_unfold:
    input: 
        T1w = bids(root='results',datatype='anat',**config['subj_wildcards'],desc='preproc', suffix='T1w.nii.gz'),
        ref = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        t1w2unfold = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='T1w',mode='image'),
    output:
        bids(root='results',datatype='anat',**config['subj_wildcards'],hemi='{hemi}',space='unfold',desc='T1w')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -d 3 -o {output} -r {input.ref} -t {input.t1w2unfold} -i {input.T1w} -v'



