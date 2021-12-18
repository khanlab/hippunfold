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
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        #warp_unfold2native_extrap = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native_extrapolateNearest.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    group: 'subj'
    log: bids(root='logs',**config['subj_wildcards'],hemi='{hemi,Lflip|R}',modality='seg-{modality}',suffix='create_warps.txt')
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
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        unfold_ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    params:
        dims = config['unfold_vol_ref']['dims'],
        script = os.path.join(config['snakemake_dir'],'workflow','scripts','mapUnfoldToFullGrid.sh'),
        warp_dir = bids(root='work',**config['subj_wildcards'],suffix='autotop',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')
    output:
        warp_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2unfoldtemplate.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfoldtemplate2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2unfoldtemplate_0InverseWarp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    container: config['singularity']['autotop']
    group: 'subj'
    threads: 8
    resources:
        time = 15 #15min
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='seg{modality}',suffix='mapUnfoldToFullGrid.txt')
    shell:
        'SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        '{params.script} {input.coords_ap} {input.coords_pd} {input.coords_io} {input.unfold_ref} {params.warp_dir} {params.dims} &> {log}'


