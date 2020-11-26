
#lookup tables for structure:
hemi_to_structure = {'L': 'HIPPOCAMPUS_LEFT', 'Lflip': 'HIPPOCAMPUS_LEFT', 'R': 'HIPPOCAMPUS_RIGHT'}
surf_to_secondary_type = {'midthickness': 'MIDTHICKNESS', 'inner': 'PIAL', 'outer': 'GRAY_WHITE'}


rule set_gifti_structure_unfold:
    input:
        unfold_gii = bids(root='work',**config['subj_wildcards'],suffix='autotop/{surfname}.unfolded.surf.gii',desc='cropped',space='{template}corobl',hemi='{hemi}',modality='{modality}'),
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname]
    output:
        unfold_gii = bids(root='work',datatype='anat',suffix='{surfname}.unfolded.surf.gii', space='{template}corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards'])
    container: config['singularity']['connectome_workbench']
    group: 'subj'
    shell:
        'wb_command -set-structure {input.unfold_gii} {params.structure_type} && ' #-surface-type FLAT -surface-secondary-type {params.secondary_type} && '
        ' cp {input} {output}'


#warp giftis to crop space
rule warp_unfold_gii_to_corobl: 
    input: 
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='{template}corobl',hemi='{hemi}',modality='{modality}'),
        unfold_gii = bids(root='work',datatype='anat',suffix='{surfname}.unfolded.surf.gii', space='{template}corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards'])
    output:
        native_gii = bids(root='work',datatype='anat',suffix='{surfname}.native.surf.gii', space='{template}corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards'])
    container: config['singularity']['connectome_workbench']
    group: 'subj'
    shell:
        'wb_command -surface-apply-warpfield {input.unfold_gii} {input.warp_unfold2native} {output.native_gii} '
#        'wb_command -set-structure {output.native_gii} 



#warp giftis to T1w space, using affine
rule warp_surf_to_T1w:
    input:
        gii = bids(root='work',datatype='anat',suffix='{surfname}.native.surf.gii', space='{template}corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='{template}corobl',desc='affineInverse',type_='ras'),
    output:
        gii = bids(root='work',datatype='anat',suffix='{surfname}.native.surf.gii', space='T1w',hemi='{hemi}',modality='{modality}', **config['subj_wildcards'],template='{template}')
    container: config['singularity']['connectome_workbench']
    group: 'subj'
    shell:
        'wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}'




