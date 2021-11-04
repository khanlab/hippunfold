import os
import numpy as np    


def get_cmd_laplace_coords(wildcards):
    if config['skip_inject_template_labels']:
        cmd = '../scripts/laplace_coords.py'
    else:
        cmd = '../scripts/laplace_coords_withinit.py'
    return cmd

def get_inputs_laplace(wildcards):
    files = dict()
    files['lbl'] = get_labels_for_laplace(wildcards)
    if not config['skip_inject_template_labels']:
        files['init_coords'] = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],dir='{dir}',suffix='coords.nii.gz',desc='init',space='corobl',hemi='{hemi}'),
    return files

rule laplace_coords:
    input: unpack(get_inputs_laplace)
    params:
        cmd = get_cmd_laplace_coords,
        gm_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['gm'],
        src_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['src'],
        sink_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['sink'],
        convergence_threshold = 1e-5,
        max_iters = 10000
    output:
        coords = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    group: 'subj'
    resources:
        time = 30
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='laplace.txt')
    script: '{params.cmd}'


rule prep_equivolume_coords:
    input: get_labels_for_laplace,
    params:
        src_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['src'],
    output:
        outerbin = bids(root='work',datatype='seg_{modality}',dir='{dir}',desc='all',suffix='mask.nii.gz',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
        innerbin = bids(root='work',datatype='seg_{modality}',dir='{dir}',desc='SRLM',suffix='mask.nii.gz',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='binarize.txt')
    group: 'subj'
    script: '../scripts/prep_equivolume_coords.py'

rule equivolume_coords:
    input: 
        outerbin = bids(root='work',datatype='seg_{modality}',dir='{dir}',desc='all',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        innerbin = bids(root='work',datatype='seg_{modality}',dir='{dir}',desc='SRLM',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    params:
        src_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['src'],
    output:
        coords = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz',desc='equivol',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    group: 'subj'
    resources:
        time = 30
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='equivolume.txt')
    container: config['singularity']['autotop'] 
    script: '../scripts/equivolume_coords.py'

rule unflip_coords:
    input:
        nii = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz', space='corobl',desc='{desc}',hemi='{hemi}flip', **config['subj_wildcards']),
    output:
        nii = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz', space='corobl',desc='{desc,laplace|equivol}',hemi='{hemi,L}', **config['subj_wildcards']),
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 'c3d {input} -flip x {output}'


#unfold ref nifti
rule create_unfold_ref:
    params:
        dims = 'x'.join(config['unfold_vol_ref']['dims']),
        voxdims = 'x'.join(config['unfold_vol_ref']['voxdims']),
        origin = 'x'.join(config['unfold_vol_ref']['origin']),
        orient = config['unfold_vol_ref']['orient']
    output: 
        nii = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}'

#this was unfold_phys_coords.nii in matlab implementation
rule create_unfold_coord_map:
    input:
        nii = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    output:
        nii = bids(root='work',space='unfold',suffix='refcoords.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input.nii} -cmp -omc {output.nii}'

def get_laminar_coords(wildcards):
        
    if 'laplace' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    elif 'equivolume' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='equivol',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    return coords_io

rule create_warps:
    input:
        unfold_ref_nii = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        unfold_phys_coords_nii = bids(root='work',space='unfold',suffix='refcoords.nii.gz',**config['subj_wildcards']),
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = get_laminar_coords
    params:
        interp_method =  'linear'
    resources:
        mem_mb = 16000
    output:
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        #warp_unfold2native_extrap = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/Warp_unfold2native_extrapolateNearest.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',modality='seg-{modality}',suffix='create_warps.txt')
    script: '../scripts/create_warps.py'


rule create_dentate_mask:
    """ threshold hipp P-D to get dentate mask
    """
    input:
        lbl = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='postproc',space='corobl',hemi='{hemi}')
    params:
        gm_labels = lambda wildcards: config['laplace_labels']['DG_IO']['gm'],
    output:
        mask = bids(root='work',datatype='seg_{modality}',label='dentate',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input.lbl} -threshold {params.gm_labels} {params.gm_labels} 1 0 {output.mask}'
        


rule laplace_dentate_io:
    input: 
        lbl = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='postproc',space='corobl',hemi='{hemi}')
    params:
        gm_labels = lambda wildcards: config['laplace_labels']['DG_IO']['gm'],
        src_labels = lambda wildcards: config['laplace_labels']['DG_IO']['src'],
        sink_labels = lambda wildcards: config['laplace_labels']['DG_IO']['sink'],
        convergence_threshold = 1e-5, 
        max_iters = 10000 
    output:
        coords = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    group: 'subj'
    resources:
        time = 30
    log: bids(root='logs',**config['subj_wildcards'],dir='IO',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='laplaceDG.txt')
    script: '../scripts/laplace_coords.py'


  
rule create_dentate_pd:
    input:
        lbl = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='postproc',space='corobl',hemi='{hemi}'),
        lbl_gm = lambda wildcards: config['laplace_labels']['DG_IO']['gm'],
        APcoords = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        APsrc_labels = lambda wildcards: config['laplace_labels'][AP]['src'],
        APsink_labels = lambda wildcards: config['laplace_labels'][AP]['sink'],
        IOcoords = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        IOsrc_labels = lambda wildcards: config['laplace_labels'][IO]['src'],
        IOsink_labels = lambda wildcards: config['laplace_labels'][IO]['sink'],
    output:
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell: '../scripts/inferGradient-CrossProd.py'


  
rule create_dentate_ap:
    input:
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        mask = bids(root='work',datatype='seg_{modality}',label='dentate',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    output:
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input.coords_ap} {input.mask} -multiply {output.coords_ap}'


       
rule create_warps_dentate:
    input: 
        unfold_ref_nii = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        unfold_phys_coords_nii = bids(root='work',space='unfold',suffix='refcoords.nii.gz',**config['subj_wildcards']),
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    params:
        interp_method =  'linear'
    resources:
        mem_mb = 16000
    output:
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopDG/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopDG/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',modality='seg-{modality}',suffix='create_warps_dentate.txt')
    script: '../scripts/create_warps.py'

       

# TODO: add this extrapolation of the warp file, extrapolate_warp_unfold2native.m
#  .. for now just using original warp
#rule extrapolate_warp_unfold2nii:
    

#full-grid correction of unfolded space
rule map_to_full_grid:
    input:
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = get_laminar_coords,
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        unfold_ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    params:
        dims = config['unfold_vol_ref']['dims'],
        script = os.path.join(config['snakemake_dir'],'workflow','scripts','mapUnfoldToFullGrid.sh'),
        warp_dir = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')
    output:
        warp_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/Warp_unfold2unfoldtemplate.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/WarpITK_unfold2unfoldtemplate_0InverseWarp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotopHipp/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    container: config['singularity']['autotop']
    group: 'subj'
    threads: 8
    resources:
        time = 15 #15min
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}',suffix='mapUnfoldToFullGrid.txt')
    shell:
        'SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        '{params.script} {input.coords_ap} {input.coords_pd} {input.coords_io} {input.unfold_ref} {params.warp_dir} {params.dims} &> {log}'

#full-grid correction of unfolded space
rule map_to_full_grid_dentate:
    input:
        coords_ap = bids(root='work',datatype='seg_{modality}',dir='AP',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_pd = bids(root='work',datatype='seg_{modality}',dir='PD',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        coords_io = bids(root='work',datatype='seg_{modality}',dir='IO',suffix='coords.nii.gz',desc='dentate',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotopDG/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        unfold_ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    params:
        dims = config['unfold_vol_ref']['dims'],
        script = os.path.join(config['snakemake_dir'],'workflow','scripts','mapUnfoldToFullGrid.sh'),
        warp_dir = bids(root='work',**config['subj_wildcards'],suffix='autotopDG',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')
    output:
        warp_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/Warp_unfold2unfoldtemplate.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/WarpITK_unfold2unfoldtemplate_0InverseWarp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotopDG/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    container: config['singularity']['autotop']
    group: 'subj'
    threads: 8
    resources:
        time = 15 #15min
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}',suffix='mapUnfoldToFullGridDG.txt')
    shell:
        'SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        '{params.script} {input.coords_ap} {input.coords_pd} {input.coords_io} {input.unfold_ref} {params.warp_dir} {params.dims} &> {log}'



rule compose_warps_corobl2unfold_rhemi:
    """ Compose corobl to unfold (unfold-template), for right hemi (ie no flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='{autotop}/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='{autotop}/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    output:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='corobl',to='unfold',mode='image',desc='{autotop}')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.unfold2unfoldtemplate} {input.native2unfold}'


rule compose_warps_corobl2unfold_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='{autotop}/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='{autotop}/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt')
    output:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='corobl',to='unfold',mode='image',desc='{autotop}')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.unfold2unfoldtemplate} {input.native2unfold} {input.flipLR_xfm}'


 
rule compose_warps_t1_to_unfold:
    """ Compose warps from T1w to unfold """
    input:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='corobl',to='unfold',mode='image',desc='{autotop}'),
        ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
    output: 
        bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='T1w',to='unfold',mode='image',desc='{autotop}')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output} -R {input.ref} {input.corobl2unfold} {input.t1w2corobl}'




