
rule import_subfield_labels:
    input: os.path.join(config['snakemake_dir'],'resources','bigbrain','sub-bigbrain_hemi-{hemi}_subfields.label.gii')
    output: bids(root='work',datatype='surf_{modality}',desc='bigbrain',suffix='subfields.labels.gii', hemi='{hemi}', **config['subj_wildcards']),
    group: 'subj'
    shell: 'cp {input} {output}'
   
 
#map bigbrain subfields to volume from surface
rule subfields_to_volume:
    input:
        label = bids(root='work',datatype='surf_{modality}',desc='bigbrain',suffix='subfields.labels.gii', hemi='{hemi}', **config['subj_wildcards']),
        midthickness = bids(root='work',datatype='surf_{modality}',suffix='midthickness.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        ref = bids(root='work',datatype='seg_{modality}',suffix='cropref.nii.gz', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        inner = bids(root='work',datatype='surf_{modality}',suffix='inner.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
        outer = bids(root='work',datatype='surf_{modality}',suffix='outer.native.surf.gii', space='T1w',hemi='{hemi}', **config['subj_wildcards']),
    output:
        vol = bids(root='work',datatype='seg_{modality}',desc='subfields',from_='surface',suffix='dseg.nii.gz', space='cropT1w',hemi='{hemi}', **config['subj_wildcards'])
    container: config['singularity']['autotop']
    group: 'subj'
    shell:
        'wb_command -label-to-volume-mapping {input.label} {input.midthickness} {input.ref} {output.vol} -ribbon-constrained {input.inner} {input.outer}'


