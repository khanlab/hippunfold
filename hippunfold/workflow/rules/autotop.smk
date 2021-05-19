import os
    

def get_input_for_shape_inject(wildcards):
    if get_modality_key(wildcards.modality) == 'seg':
        modality_suffix = get_modality_suffix(wildcards.modality)
        seg = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',from_='{modality_suffix}').format(
                    **wildcards, modality_suffix=modality_suffix),
    else:
        seg = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='nnunet',space='corobl',hemi='{hemi}').format(**wildcards)
    return seg

rule template_shape_inject:
    """should refactor into snakemake so it doesn't require fsl+ants in container"""
    input:
        template_seg = os.path.join(config['snakemake_dir'],'resources','UPenn_templateShape','labelmap.nii.gz'),
        seg = get_input_for_shape_inject
    params:
        preserve_labels = 7,
        script = os.path.join(config['snakemake_dir'],'workflow','scripts','template_shape_inject.sh')
    output:
        postproc_tmpdir = directory(bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}')),
        postproc = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess.nii.gz',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
    threads: 8
    resources:
        time = 180 #3 hrs (so grouped job is set at 3hrs, since snakemake grouped resources not summed, but takes the max)
    group: 'subj'
    container: config['singularity']['autotop'] #currently requires fsl and ants
    log: bids(root='logs',**config['subj_wildcards'],space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}',suffix='templateshapereg.txt')
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        '{params.script} {input.template_seg} {input.seg} {output.postproc} {output.postproc_tmpdir} -L {params.preserve_labels} &> {log}' 



rule transform_init_coords:
    """ not currently used"""
    input:
        postproc_tmpdir = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess', desc='cropped', space='corobl', hemi='{hemi}', modality='{modality}'),
        ref = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess.nii.gz', desc='cropped', space='corobl', hemi='{hemi}', modality='{modality}'),
        coords = os.path.join(config['snakemake_dir'],'resources','UPenn_templateShape','coords-{dir}.nii.gz')
    params:
        in_aff = lambda wildcards, input: os.path.join(input.postproc_tmpdir,'ants_0GenericAffine.mat'),
        in_warp = lambda wildcards, input: os.path.join(input.postproc_tmpdir,'ants_1Warp.nii.gz'),
        interp = 'NearestNeighbor'
    output:
        coords = bids(root='work',suffix='autotop/initcoords-{dir}.nii.gz',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards'])
    group: 'subj'
    container: config['singularity']['ants']
    shell:
        'ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS={threads} '
        "antsApplyTransforms -d 3 -i {input.coords} -r {input.ref} -o {output.coords} -n {params.interp} -t {params.in_warp} -t {params.in_aff} {params.autotop_cmd}"

rule laplace_coords:
    input:
        lbl = bids(root='work',**config['subj_wildcards'],suffix='autotop/labelmap-postProcess.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    output:
        coords_ap = bids(root='work',suffix='autotop/coords-AP.nii.gz',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),
        coords_pd = bids(root='work',suffix='autotop/coords-PD.nii.gz',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),
        coords_io = bids(root='work',suffix='autotop/coords-IO.nii.gz',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),
    group: 'subj'
    script: '../scripts/laplace_coords.py'


rule placeholder_warps_gifti:
    input:
        coords_ap = bids(root='work',suffix='autotop/coords-AP.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        coords_pd = bids(root='work',suffix='autotop/coords-PD.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        coords_io = bids(root='work',suffix='autotop/coords-IO.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
    output:
        warp_unfold2native_extrap = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native_extrapolateNearest.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warp_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        warpitk_native2unfold= bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}'),
        gii = expand(bids(root='work',suffix='autotop/{surfname}.unfoldedtemplate.surf.gii',desc='cropped', space='corobl',hemi='{hemi,Lflip|R}',modality='{modality}', **config['subj_wildcards']),surfname=['inner','outer','midthickness'],allow_missing=True),
    group: 'subj'
    shell: 'echo placeholder only'


#full-grid correction of unfolded space
rule map_to_full_grid:
    input:
        coords_ap = bids(root='work',suffix='autotop/coords-AP.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        coords_pd = bids(root='work',suffix='autotop/coords-PD.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        coords_io = bids(root='work',suffix='autotop/coords-IO.nii.gz',desc='cropped', space='corobl',hemi='{hemi}',modality='{modality}', **config['subj_wildcards']),
        warp_unfold2native = bids(root='work',**config['subj_wildcards'],suffix='autotop/Warp_unfold2native.nii',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
        unfold_ref = os.path.join(config['snakemake_dir'],'hippocampal_autotop','misc','unfold_ref_256x128x16.nii.gz'),
    params:
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
        '{params.script} {input.coords_ap} {input.coords_pd} {input.coords_io} {input.unfold_ref} {params.warp_dir} &> {log}'


#rule to unflip a nifti
rule unflip_autotop_nii:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/{filename}.nii.gz',desc='cropped',space='corobl',hemi='{hemi}flip',modality='{modality}')
    output:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/{filename}.nii.gz',desc='cropped',space='corobl',hemi='{hemi,L}',modality='{modality}')
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 'c3d {input} -flip x {output}'


rule cp_coords_to_bids:
    input:
        nii = bids(root='work',**config['subj_wildcards'],suffix='autotop/coords-{dir}.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',modality='{modality}'),
    output:
        nii = bids(root='work',datatype='seg_{modality}',dir='{dir}',suffix='coords.nii.gz', space='corobl',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 'cp {input} {output}'



rule compose_warps_corobl2unfold_rhemi:
    """ Compose corobl to unfold (unfold-template), for right hemi (ie no flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='R',modality='{modality}'),
        ref = os.path.join(config['snakemake_dir'],'hippocampal_autotop','misc','unfold_ref_256x128x16.nii.gz'),
    output:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='R',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.unfold2unfoldtemplate} {input.native2unfold}'


rule compose_warps_corobl2unfold_lhemi:
    """ Compose corobl to unfold (unfold-template), for left hemi (ie with flip)"""
    input:
        native2unfold = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_native2unfold.nii',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        unfold2unfoldtemplate = bids(root='work',**config['subj_wildcards'],suffix='autotop/WarpITK_unfold2unfoldtemplate_0Warp.nii.gz',desc='cropped',space='corobl',hemi='Lflip',modality='{modality}'),
        ref = os.path.join(config['snakemake_dir'],'hippocampal_autotop','misc','unfold_ref_256x128x16.nii.gz'),
        flipLR_xfm = os.path.join(config['snakemake_dir'],'resources','desc-flipLR_type-itk_xfm.txt')
    output:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='L',from_='corobl',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output.corobl2unfold} -R {input.ref} {input.unfold2unfoldtemplate} {input.native2unfold} {input.flipLR_xfm}'


 
rule compose_warps_t1_to_unfold:
    """ Compose warps from T1w to unfold """
    input:
        corobl2unfold = bids(root='work',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='corobl',to='unfold',mode='image'),
        ref = os.path.join(config['snakemake_dir'],'hippocampal_autotop','misc','unfold_ref_256x128x16.nii.gz'),
        t1w2corobl = bids(root='work',datatype='anat',**config['subj_wildcards'],suffix='xfm.txt',from_='T1w',to='corobl',desc='affine',type_='itk'),
    output: 
        bids(root='results',datatype='seg_{modality}',**config['subj_wildcards'],suffix='xfm.nii.gz',hemi='{hemi}',from_='T1w',to='unfold',mode='image')
    container: config['singularity']['ants']
    group: 'subj'
    shell: 'ComposeMultiTransform 3 {output} -R {input.ref} {input.corobl2unfold} {input.t1w2corobl}'



