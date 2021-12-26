rule create_native_coord_ref:
    input:
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    output:
        nii = bids(root='work',datatype='seg_{modality}',space='corobl',hemi='{hemi}',suffix='refcoords.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input} -cmp -omc {output}'



rule create_warps:
    input:
        unfold_ref_nii = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refvol.nii.gz',**config['subj_wildcards']),
        unfold_phys_coords_nii = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refcoords.nii.gz',**config['subj_wildcards']),
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = get_laminar_coords,
        native_ref_coords_nii = bids(root='work',datatype='seg_{modality}',space='corobl',hemi='{hemi}',suffix='refcoords.nii.gz',**config['subj_wildcards'])
    params:
        interp_method =  'linear'
    resources:
        mem_mb = 16000
    output:
        warp_unfold2native = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='surface'),
        warp_native2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='surface'),
        warpitk_unfold2native = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='image'),

        warpitk_native2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='image')
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',modality='seg-{modality}',suffix='create_warps.txt')
    script: '../scripts/create_warps.py'



rule compose_warps_corobl2unfold_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='Lflip',from_='corobl',to='unfold',mode='image'),
        ref = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refvol.nii.gz',**config['subj_wildcards']),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt')
    output:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.native2unfold} {input.flipLR_xfm}'


#consider renaming to state this composition is to "un-flip"
rule compose_warps_unfold2corobl_lhemi:
    """ Compose unfold (unfold-template) to corobl, for left hemi (ie with flip)"""
    input:
        unfold2native = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='Lflip',from_='unfold',to='corobl',mode='image'),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt'),
        unfold_ref = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refvol.nii.gz',**config['subj_wildcards']),
        ref = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='L', **config['subj_wildcards']),
    output:
        unfold2corobl = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='unfold',to='corobl',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2corobl},1] -r {input.ref} -t {input.flipLR_xfm} -t {input.unfold2native} -i {input.unfold_ref} -v'


 
rule compose_warps_t1_to_unfold:
    """ Compose warps from T1w to unfold """
    input:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='corobl',to='unfold',mode='image'),
        ref = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refvol.nii.gz',**config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
    output: 
        bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='T1w',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output} -R {input.ref} {input.corobl2unfold} {input.t1w2corobl}'


rule compose_warps_unfold_to_cropt1:
    """ Compose warps from unfold to cropT1w """
    input:
        unfold2corobl = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='corobl',mode='image'),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        unfold_ref = bids(root='work',space='unfold',datatype='seg_{modality}',suffix='refvol.nii.gz',**config['subj_wildcards']),
    output: 
        unfold2cropt1 = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='unfold',to='T1w',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2cropt1},1] -r {input.ref} -t [{input.t1w2corobl},1] -t {input.unfold2corobl} -i {input.unfold_ref} -v'




