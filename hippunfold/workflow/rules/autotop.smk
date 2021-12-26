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
        files['init_coords'] = bids(root='work',datatype='seg',**config['subj_wildcards'],dir='{dir}',suffix='coords.nii.gz',desc='init',space='corobl',hemi='{hemi}'),
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
        coords = bids(root='work',datatype='seg',dir='{dir}',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    group: 'subj'
    resources:
        time = 30
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',suffix='laplace.txt')
    script: '{params.cmd}'

rule prep_equivolume_coords:
    input: get_labels_for_laplace,
    params:
        src_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['src'],
    output:
        outerbin = bids(root='work',datatype='seg',dir='{dir}',desc='all',suffix='mask.nii.gz',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
        innerbin = bids(root='work',datatype='seg',dir='{dir}',desc='SRLM',suffix='mask.nii.gz',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',suffix='binarize.txt')
    group: 'subj'
    script: '../scripts/prep_equivolume_coords.py'

rule equivolume_coords:
    input: 
        outerbin = bids(root='work',datatype='seg',dir='{dir}',desc='all',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        innerbin = bids(root='work',datatype='seg',dir='{dir}',desc='SRLM',suffix='mask.nii.gz',space='corobl',hemi='{hemi}', **config['subj_wildcards']),
    params:
        src_labels = lambda wildcards: config['laplace_labels'][wildcards.dir]['src'],
    output:
        coords = bids(root='work',datatype='seg',dir='{dir}',suffix='coords.nii.gz',desc='equivol',space='corobl',hemi='{hemi,Lflip|R}', **config['subj_wildcards']),
    group: 'subj'
    resources:
        time = 30
    log: bids(root='logs',**config['subj_wildcards'],dir='{dir}',hemi='{hemi,Lflip|R}',suffix='equivolume.txt')
    container: config['singularity']['autotop'] 
    script: '../scripts/equivolume_coords.py'

rule unflip_coords:
    input:
        nii = bids(root='work',datatype='seg',dir='{dir}',suffix='coords.nii.gz', space='corobl',desc='{desc}',hemi='{hemi}flip', **config['subj_wildcards']),
    output:
        nii = bids(root='work',datatype='seg',dir='{dir}',suffix='coords.nii.gz', space='corobl',desc='{desc,laplace|equivol}',hemi='{hemi,L}', **config['subj_wildcards']),
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
        nii = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d -create {params.dims} {params.voxdims}mm -origin {params.origin}mm -orient {params.orient} -o {output.nii}'

#this was unfold_phys_coords.nii in matlab implementation
rule create_unfold_coord_map:
    input:
        nii = bids(root='work',datatype='seg',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards'])
    output:
        nii = bids(root='work',datatype='seg',space='unfold',suffix='refcoords.nii.gz',**config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop'] 
    shell:
        'c3d {input.nii} -cmp -omc {output.nii}'

def get_laminar_coords(wildcards):
        
    if 'laplace' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg',dir='IO',suffix='coords.nii.gz',desc='laplace',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    elif 'equivolume' in config['laminar_coords_method']:
        coords_io = bids(root='work',datatype='seg',dir='IO',suffix='coords.nii.gz',desc='equivol',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    return coords_io


