
rule t1_t2_ratio:
    input:
       t1=bids(
            root=root,
            datatype=anat,
            space='{space}',
            desc='preproc',
            suffix='T1w',
            **config["subj_wildcards"])
        t2=bids(
            root=root,
            datatype=anat,
            space='{space}',
            desc='preproc',
            suffix='T2w',
            **config["subj_wildcards"])
    output:
        t1overt2=bids(
            root=root,
            datatype=anat,
            space='{space}',
            desc='preproc',
            suffix='T1wT2wRatio',
            **config["subj_wildcards"])
    group: 'subj'
    container:
        config["singularity"]["autotop"]
    shell:
        'c3d {input.t1} {input.t2} -divide -replace inf 1000 -inf -1000 NaN 0' 

#sample on hipp & dg midthickness surfaces

                 
