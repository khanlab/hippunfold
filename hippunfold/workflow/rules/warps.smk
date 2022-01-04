rule create_native_coord_ref:
    input:
        coords_ap = bids(root='work',datatype='seg',dir='AP',suffix='coords-{autotop}.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    output:
        nii = bids(root='work',datatype='seg',space='corobl',hemi='{hemi}',suffix='refcoords-{autotop}.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input} -cmp -omc {output}'

#unfold ref nifti
rule create_unfold_ref:
    params:
        dims = lambda wildcards: 'x'.join(config[wildcards.autotop]['dims']),
        voxdims = lambda wildcards: 'x'.join(config[wildcards.autotop]['voxdims']),
        origin = lambda wildcards: 'x'.join(config[wildcards.autotop]['origin']),
        orient = lambda wildcards: config[wildcards.autotop]['orient']
    output: 
        nii = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}'

#this was unfold_phys_coords.nii in matlab implementation
rule create_unfold_coord_map:
    input:
        nii = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards'])
    output:
        nii = bids(root='work',space='unfold',suffix='refcoords-{autotop}.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input.nii} -cmp -omc {output.nii}'


def get_laminar_coords(wildcards):
    if 'laplace' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg',dir='IO',suffix='coords-autotopHipp.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    elif 'equivolume' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg',dir='IO',suffix='coords-autotopHipp.nii.gz',desc='equivol',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    return coords_io

rule create_warps_hipp:
    input:
        unfold_ref_nii = bids(root='work',space='unfold',suffix='refvol-autotopHipp.nii.gz',**config['subj_wildcards']),
        unfold_phys_coords_nii = bids(root='work',space='unfold',suffix='refcoords-autotopHipp.nii.gz',**config['subj_wildcards']),
        coords_ap = bids(root='work',datatype='seg',dir='AP',suffix='coords-autotopHipp.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg',dir='PD',suffix='coords-autotopHipp.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = get_laminar_coords,
        native_ref_coords_nii = bids(root='work',datatype='seg',space='corobl',hemi='{hemi}',suffix='refcoords-autotopHipp.nii.gz',**config['subj_wildcards'])
    params:
        interp_method =  'linear'
    resources:
        mem_mb = 16000
    output:
        warp_unfold2native = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopHipp.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='surface'),
        warp_native2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopHipp.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='surface'),
        warpitk_unfold2native = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopHipp.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='image'),
        warpitk_native2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopHipp.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='image')
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',suffix='create_warps-autotopHipp.txt') 
    script: '../scripts/create_warps.py'


rule create_warps_DG:
    input:
        unfold_ref_nii = bids(root='work',space='unfold',suffix='refvol-autotopDG.nii.gz',**config['subj_wildcards']),
        unfold_phys_coords_nii = bids(root='work',space='unfold',suffix='refcoords-autotopDG.nii.gz',**config['subj_wildcards']),
        coords_ap = bids(root='work',datatype='seg',dir='AP',suffix='coords-autotopDG.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg',dir='PD',suffix='coords-autotopDG.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = bids(root='work',datatype='seg',dir='IO',suffix='coords-autotopDG.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        native_ref_coords_nii = bids(root='work',datatype='seg',space='corobl',hemi='{hemi}',suffix='refcoords-autotopDG.nii.gz',**config['subj_wildcards'])
    params:
        interp_method =  'linear'
    resources:
        mem_mb = 16000
    output:
        warp_unfold2native = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopDG.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='surface'),
        warp_native2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopDG.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='surface'),
        warpitk_unfold2native = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopDG.nii.gz',hemi='{hemi,Lflip|R}',from_='unfold',to='corobl',mode='image'),
        warpitk_native2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-autotopDG.nii.gz',hemi='{hemi,Lflip|R}',from_='corobl',to='unfold',mode='image')
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',suffix='create_warps-autotopDG.txt') 
    script: '../scripts/create_warps.py'



rule compose_warps_corobl2unfold_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='Lflip',from_='corobl',to='unfold',mode='image'),
        ref = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards']),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt')
    output:
        corobl2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='L',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.native2unfold} {input.flipLR_xfm}'


#consider renaming to state this composition is to "un-flip"
rule compose_warps_unfold2corobl_lhemi:
    """ Compose unfold (unfold-template) to corobl, for left hemi (ie with flip)"""
    input:
        unfold2native = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='Lflip',from_='unfold',to='corobl',mode='image'),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt'),
        unfold_ref = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards']),
        ref = bids(root='work',datatype='seg',dir='AP',suffix='coords-{autotop}.nii.gz',desc='laplace',space='corobl',hemi='L', **config['subj_wildcards']),
    output:
        unfold2corobl = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='L',from_='unfold',to='corobl',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2corobl},1] -r {input.ref} -t {input.flipLR_xfm} -t {input.unfold2native} -i {input.unfold_ref} -v'


 
rule compose_warps_t1_to_unfold:
    """ Compose warps from T1w to unfold """
    input:
        corobl2unfold = bids(root='work',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='{hemi}',from_='corobl',to='unfold',mode='image'),
        ref = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
    output: 
        bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='{hemi}',from_='T1w',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output} -R {input.ref} {input.corobl2unfold} {input.t1w2corobl}'


rule compose_warps_unfold_to_cropt1:
    """ Compose warps from unfold to cropT1w """
    input:
        unfold2corobl = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='{hemi}',from_='unfold',to='corobl',mode='image'),
        ref = bids(root='work',datatype='seg',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        unfold_ref = bids(root='work',space='unfold',suffix='refvol-{autotop}.nii.gz',**config['subj_wildcards']),
    output: 
        unfold2cropt1 = bids(root='results',datatype='seg',**config['subj_wildcards'],suffix='xfm-{autotop}.nii.gz',hemi='{hemi}',from_='unfold',to='T1w',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'antsApplyTransforms -o [{output.unfold2cropt1},1] -r {input.ref} -t [{input.t1w2corobl},1] -t {input.unfold2corobl} -i {input.unfold_ref} -v'




