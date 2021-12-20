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



rule compose_warps_rhemi:
    """ Compose corobl to unfold (unfold-template), for right hemi (ie no flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        refnative = bids(root='results',datatype='anat',**config['subj_wildcards'],desc='preproc',suffix='T1w.nii.gz')
    output:
        unfold2native = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='T1w',to='unfold',mode='image'),
        native2unfold = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='unfold',to='T1w',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.unfold2native} -R {input.ref} {input.t1w2corobl} {input.native2unfold} && '
            'ComposeMultiTransform 3 {output.native2unfold} -R {input.refnative} -i {input.t1w2corobl} {input.unfold2native}'


rule compose_warps_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
        ref = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt'),
        refnative = bids(root='results',datatype='anat',**config['subj_wildcards'],desc='preproc',suffix='T1w.nii.gz')
    output:
        unfold2native = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='T1w',to='unfold',mode='image'),
        native2unfold = bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='unfold',to='T1w',mode='image')
    container: config['singularity']['autotop'] 
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.unfold2native} -R {input.ref} {input.t1w2corobl} {input.native2unfold} {input.flipLR_xfm} && '
            'ComposeMultiTransform 3 {output.native2unfold} -R {input.refnative} -i {input.t1w2corobl} {input.unfold2native} -i {input.flipLR_xfm}'


