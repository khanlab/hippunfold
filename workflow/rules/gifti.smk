
#lookup tables for structure:
hemi_to_structure = {'L': 'CORTEX_LEFT', 'Lflip': 'CORTEX_LEFT', 'R': 'CORTEX_RIGHT'}
surf_to_secondary_type = {'midthickness': 'MIDTHICKNESS', 'inner': 'PIAL', 'outer': 'GRAY_WHITE'}


rule warp_gii_unfoldtemplate2unfold: 
    """warp from template space to subj unfolded"""
    input: 
        warp = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfoldtemplate2unfold.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        gii = os.path.join(config['snakemake_dir'],'resources','unfold_template','tpl-avg_space-unfold_den-{density}_{surfname}.surf.gii')
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type = 'FLAT'
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded', hemi='{hemi,R|Lflip}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && '
        'wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}'
            ' -surface-secondary-type {params.secondary_type}'


rule calc_unfold_template_coords:
    """ Creates a coords.shape.gii from the unfolded template """
    input:
        midthickness_gii = os.path.join(config['snakemake_dir'],'resources','unfold_template','tpl-avg_space-unfold_den-{density}_midthickness.surf.gii')
    params:
        coords_xyz = 'coords-XYZ.shape.gii',
        coord_AP = 'coord-AP.shape.gii',
        coord_PD = 'coord-PD.shape.gii',
        coord_IO = 'coord-IO.shape.gii',
        origin = config['unfold_vol_ref']['origin'],
        extent = config['unfold_vol_ref']['extent'],
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
    output:
        coords_gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='coords.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
    container: config['singularity']['autotop']
    shadow: 'minimal' #this is required to use the temporary files defined as params
    group: 'subj'
    shell:
        "wb_command -surface-coordinates-to-metric {input.midthickness_gii} {params.coords_xyz} && "
        "wb_command -file-information {params.coords_xyz} && "
        "wb_command -metric-math '({params.origin[0]}-X)/{params.extent[0]}' {params.coord_AP} -var X {params.coords_xyz} -column 1 && "
        "wb_command -file-information {params.coord_AP} && "
        "wb_command -metric-math '({params.origin[1]}+Y)/{params.extent[1]}' {params.coord_PD} -var Y {params.coords_xyz} -column 2 && "
        "wb_command -file-information {params.coord_PD} && "
        "wb_command -metric-math '({params.origin[2]}+Z)/{params.extent[2]}' {params.coord_IO} -var Z {params.coords_xyz} -column 3 && "
        "wb_command -file-information {params.coord_IO} && "
        "wb_command -metric-merge {output.coords_gii}  -metric {params.coord_AP} -metric {params.coord_PD} -metric {params.coord_IO} && "
        "wb_command -set-structure {output.coords_gii} {params.structure_type}"

       
        
   

    

#subj unfolded surf might have a few vertices outside the bounding box.. this constrains all the vertices to the warp bounding box
rule constrain_surf_to_bbox:
    input:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded',hemi='{hemi}', **config['subj_wildcards']),
        ref_nii = bids(root='work',space='unfold',suffix='refvol.nii.gz',**config['subj_wildcards']),
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii',desc='constrainbbox', space='unfolded',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    script: '../scripts/constrain_surf_to_bbox.py'

#warp from subj unfolded to corobl
rule warp_gii_unfold2native: 
    input: 
        warp = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', desc='constrainbbox',space='unfolded',hemi='{hemi}', **config['subj_wildcards'])
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type = 'ANATOMICAL'
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', desc='nonancorrect', space='corobl',hemi='{hemi,R|Lflip}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-warpfield {input.gii} {input.warp} {output.gii} && '
        'wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}'
            ' -surface-secondary-type {params.secondary_type}'


# previous rule seems to be where nan vertices emerge, so we'll correct them here immediately after
rule correct_nan_vertices:
    input: 
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', desc='nonancorrect', space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='corobl',hemi='{hemi,R|Lflip}', **config['subj_wildcards'])
    group: 'subj'
    script: '../scripts/fillnanvertices.py'


#unflip surface
rule unflip_gii:
    input:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='corobl',hemi='{hemi}flip', **config['subj_wildcards'])
    params:
        structure_type = lambda wildcards: hemi_to_structure[wildcards.hemi],
        secondary_type = lambda wildcards: surf_to_secondary_type[wildcards.surfname],
        surface_type = 'ANATOMICAL'
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='corobl',hemi='{hemi,L}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-flip-lr {input.gii} {output.gii} && '
        'wb_command -set-structure {output.gii} {params.structure_type} -surface-type {params.surface_type}'
        ' -surface-secondary-type {params.secondary_type}'



def get_unfolded_surf_R_Lflip (wildcards):
    if wildcards.hemi == 'R':
        return  bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded',hemi='{hemi}', **config['subj_wildcards']).format(**wildcards)
    elif wildcards.hemi == 'L':
        return  bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded',hemi='{hemi}flip', **config['subj_wildcards']).format(**wildcards)


rule unflip_gii_unfolded:
    """copy unfolded from Lflip to L"""
    input:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded',hemi='{hemi}flip', **config['subj_wildcards'])
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='unfolded',hemi='{hemi,L}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'cp {input.gii} {output.gii}'



#warp from corobl to T1w
rule warp_gii_to_T1w:
    input:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='corobl',hemi='{hemi}', **config['subj_wildcards']),
        xfm = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='ras'),
    output:
        gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -surface-apply-affine {input.gii} {input.xfm} {output.gii}'



#morphological features, calculated in T1w space:
rule calculate_surface_area:
    input: 
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='midthickness.surf.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='surfarea.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-vertex-areas {input} {output}"    

rule calculate_gyrification:
    """new gyrification is ratio of nativearea to unfoldarea (e.g. surface scaling or distortion factor.
    this should be proportional by a constant, to the earlier gyrification on 32k surfaces."""
    input: 
        native_surfarea = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='surfarea.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
        unfold_surfarea = os.path.join(config['snakemake_dir'],'resources','unfold_template','tpl-avg_space-unfold_den-{density}_surfarea.shape.gii')
    output:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='gyrification.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        'wb_command -metric-math "nativearea/unfoldarea" {output.gii}'
            ' -var nativearea {input.native_surfarea} -var unfoldarea {input.unfold_surfarea}'


rule smooth_surface:
    input:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='midthickness.surf.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    params:
        strength = 0.6,
        iterations = 100
    output:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='midthickness.surf.gii', space='{space}',hemi='{hemi}', desc='smoothed',**config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-smoothing {input.gii} {params.strength} {params.iterations} {output.gii}"

    

rule calculate_curvature_from_surface:
    input: 
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='midthickness.surf.gii', space='{space}',hemi='{hemi}', desc='smoothed',**config['subj_wildcards'])
    output:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='curvature.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-curvature {input} -mean {output}"    

        
rule calculate_thickness_from_surface:
    input: 
        inner = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='inner.surf.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
        outer = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='outer.surf.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    output:
        gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='thickness.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell:
        "wb_command -surface-to-surface-3d-distance {input.outer} {input.inner} {output}"    

rule resample_bigbrain_subfield_label_gii:
    """ similar to wb_command -metric-resample, but for unfolded space
         note, this creates a Nx1 nii file, for subsequent gifti conversion
         using wb_command""" 
    input:
        label = os.path.join(config['snakemake_dir'],'resources','bigbrain','sub-bigbrain_hemi-{hemi}_subfields.label.gii'),
        new_surf = os.path.join(config['snakemake_dir'],'resources','unfold_template','tpl-avg_space-unfold_den-{density}_midthickness.surf.gii')
    output:
        label_nii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='subfields.label.nii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    script: '../scripts/resample_unfolded_label.py'

rule label_nii_to_metric_gii:
    input:
        label_nii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='subfields.label.nii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
        surf = os.path.join(config['snakemake_dir'],'resources','unfold_template','tpl-avg_space-unfold_den-{density}_midthickness.surf.gii')
    output:
        metric_gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='subfields.label.func.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop']
    shell: 'wb_command -metric-convert -from-nifti {input.label_nii} {input.surf} {output.metric_gii}'

rule metric_to_label_gii:
    input:
        metric_gii = bids(root='work',datatype='surf_{modality}',den='{density}',suffix='subfields.label.func.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
        label_list = os.path.join(config['snakemake_dir'],'resources','bigbrain','sub-bigbrain_labellist.txt')
    output:
        label_gii = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='subfields.label.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['autotop']
    shell: 'wb_command -metric-label-import {input.metric_gii} {input.label_list} {output.label_gii}'


def get_cmd_cifti_metric(wildcards,input, output):
    cmd = f'wb_command  -cifti-create-dense-scalar {output}'
    if 'L' in config['hemi']:
        cmd = cmd + f' -left-metric {input.left_metric}'
    if 'R' in config['hemi']:
        cmd = cmd + f' -right-metric {input.right_metric}'
    return cmd

def get_inputs_cifti_metric(wildcards):
    files = dict()
    if 'L' in config['hemi']:
        files['left_metric'] = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{metric}.shape.gii', space='{space}',hemi='L', **config['subj_wildcards']).format(**wildcards),
    if 'R' in config['hemi']:
        files['right_metric'] = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{metric}.shape.gii', space='{space}',hemi='R', **config['subj_wildcards']).format(**wildcards),
    return files

rule create_dscalar_metric_cifti:
    input: unpack(get_inputs_cifti_metric)
    params:
        cmd = get_cmd_cifti_metric
    output:
        cifti = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{metric}.dscalar.nii', space='{space}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell: '{params.cmd}'

def get_inputs_cifti_label(wildcards):
    files = dict()
    if 'L' in config['hemi']:
        files['left_label'] = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='subfields.label.gii', space='{space}',hemi='L', **config['subj_wildcards']).format(**wildcards),
    if 'R' in config['hemi']:
        files['right_label'] = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='subfields.label.gii', space='{space}',hemi='R', **config['subj_wildcards']).format(**wildcards),
    return files


def get_cmd_cifti_label(wildcards,input, output):
    cmd = f'wb_command  -cifti-create-label {output}'
    if 'L' in config['hemi']:
        cmd = cmd + f' -left-label {input.left_label}'
    if 'R' in config['hemi']:
        cmd = cmd + f' -right-label {input.right_label}'
    return cmd


rule create_dlabel_cifti_subfields:
    input: unpack(get_inputs_cifti_label)
    params:
        cmd = get_cmd_cifti_label
    output:
        cifti = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='subfields.dlabel.nii', space='{space}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell: '{params.cmd}'



def get_cmd_spec_file(wildcards, input, output):
    specfile = output.spec_file
    if 'hemi' in wildcards._names:
        structure = hemi_to_structure[wildcards.hemi]
    else:
        structure = 'INVALID'
    cmds = list()
    for infile in input:
        cmds.append(' '.join(['wb_command','-add-to-spec-file',specfile, structure, infile]))
    return ' && '.join(cmds)

    

#add surfs and metrics to a spec file
rule create_spec_file:
    input:
        shapes = expand(bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{shape}.shape.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
                    shape=['gyrification','curvature','thickness','coords'], allow_missing=True),
        surfs = expand(bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{surfname}.surf.gii', space='{space}', hemi='{hemi}', **config['subj_wildcards']),
                    surfname=['midthickness','inner','outer'], space=['{space}','unfolded'], allow_missing=True), 
        subfields = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='subfields.label.gii', space='{space}',hemi='{hemi}', **config['subj_wildcards']),
        cifti = expand(bids(root='results',datatype='surf_{modality}',den='{density}',suffix='{cifti}.nii', space='{space}', **config['subj_wildcards']),
                    cifti=['gyrification.dscalar','curvature.dscalar','thickness.dscalar','subfields.dlabel'], allow_missing=True),
 
    params:
        cmds = get_cmd_spec_file
    output: 
        spec_file = bids(root='results',datatype='surf_{modality}',den='{density}',suffix='hippunfold.spec', hemi='{hemi,L|R}',space='{space}',**config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj' 
    shell: '{params.cmds}'



rule merge_lr_spec_file:
    input:
        spec_files = expand(bids(root='results',datatype='surf_{modality}',den='{density}',suffix='hippunfold.spec', hemi='{hemi}',space='{space}',**config['subj_wildcards']),
                        hemi=['L','R'], allow_missing=True)
    output:
        spec_file = bids(root='work',datatype='surf_{modality}',den='{density}',space='{space}',suffix='hippunfold.spec', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 'wb_command -spec-file-merge {input.spec_files} {output}'





