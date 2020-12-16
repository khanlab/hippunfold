
#lookup tables for structure:
hemi_to_structure = {'L': 'CORTEX_LEFT', 'Lflip': 'CORTEX_LEFT', 'R': 'CORTEX_RIGHT'}
surf_to_secondary_type = {'midthickness': 'MIDTHICKNESS', 'inner': 'PIAL', 'outer': 'GRAY_WHITE'}




#warp from template space to subj unfolded
rule warp_gii_unfoldtemplate2unfold: 
    input: 
        warp = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        gii = bids(root='work',suffix='autotop/{surfname}.unfoldedtemplate.surf.gii',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards'])
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type = 'FLAT'
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', space='corobl', hemi='{hemi,R|Lflip}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && '
        'wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}'
            ' -surface-secondary-type {params.secondary_type}'

#subj unfolded surf might have a few vertices outside the bounding box.. this constrains all the vertices to the warp bounding box
rule constrain_surf_to_bbox:
    input:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        ref_nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/unfold_ref_256x128x16.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii',desc='constrainbbox', space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    script: '../scripts/constrain_surf_to_bbox.py'

#warp from subj unfolded to corobl
rule warp_gii_unfold2native: 
    input: 
        warp = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native_extrapolateNearest.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', desc='constrainbbox',space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type = 'ANATOMICAL'
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='corobl',hemi='{hemi,R|Lflip}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && '
        'wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}'
            ' -surface-secondary-type {params.secondary_type}'

#unflip surface
rule unflip_gii:
    input:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='corobl',hemi='{hemi}flip', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='corobl',hemi='{hemi,L}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-flip-lr {input.gii} {output.gii}'

#for Lflip to L unfolded, we just copy as is without moving vertices: 
rule cp_to_unflip_unfolded:
    input:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', space='corobl',hemi='{hemi}flip', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', space='corobl',hemi='{hemi,L}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'cp {input.gii} {output.gii}'




#warp from corobl to T1w
rule warp_gii_to_T1w:
    input:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='ras'),
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}'


#morphological features, calculated in T1w space:
rule calculate_gyrification_from_surface:
    input: 
        gii = bids(root='work',datatype='surf_{modality}',suffix='midthickness.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='gyrification.native.shape.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    


rule calculate_curvature_from_surface:
    input: 
        gii = bids(root='work',datatype='surf_{modality}',suffix='midthickness.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='curvature.native.shape.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-curvature {input} -mean {output}"    

        
rule calculate_thickness_from_surface:
    input: 
        inner = bids(root='work',datatype='surf_{modality}',suffix='inner.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        outer = bids(root='work',datatype='surf_{modality}',suffix='outer.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',suffix='thickness.native.shape.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"    

def get_cmd_spec_file(wildcards, input, output):
    specfile = output.spec_file
    structure = hemi_to_structure[wildcards.hemi]
    cmds = list()
    for infile in input:
        cmds.append(' '.join(['wb_command','-add-to-spec-file',specfile, structure, infile]))
    return ' && '.join(cmds)

    

#add surfs and metrics to a spec file
rule create_spec_file:
    input:
        shapes = expand(bids(root='work',datatype='surf_{modality}',suffix='{shape}.native.shape.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
                    shape=['gyrification','curvature','thickness'], allow_missing=True),
        surfs = expand(bids(root='work',datatype='surf_{modality}',suffix='{surfname}.native.surf.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
                    surfname=['midthickness','inner','outer'], space=['T1w','corobl'], allow_missing=True), 
        unfolded = expand(bids(root='work',datatype='surf_{modality}',suffix='{surfname}.unfolded.surf.gii', space='corobl',hemi='{hemi}', **config['subj_wildcards']),
                    surfname=['midthickness','inner','outer'],  allow_missing=True), 
    params:
        cmds = get_cmd_spec_file
    output: 
        spec_file = bids(root='work',datatype='surf_{modality}',suffix='hippunfold.spec', hemi='{hemi,L|R}',**config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell: '{params.cmds}'

rule merge_lr_spec_file:
    input:
        spec_files = expand(bids(root='work',datatype='surf_{modality}',suffix='hippunfold.spec', hemi='{hemi}',**config['subj_wildcards']),
                        hemi=['L','R'], allow_missing=True)
    output: 
        spec_file = bids(root='work',datatype='surf_{modality}',suffix='hippunfold.spec', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell: 'wb_command -spec-file-merge {input.spec_files} {output}'


