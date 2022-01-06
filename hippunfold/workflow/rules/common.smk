
#take mean of all scans if >1, otherwise just copy the one scan
def get_avg_or_cp_scans_cmd (wildcards, input, output):
    if len(input) > 1:
        cmd = f'c3d {input} -mean -o {output}'
    else:
        cmd = f'cp {input} {output}'
    return cmd



def get_modality_key(modality):
    if modality[:3] == 'seg': 
        return 'seg'
    else:
        return modality

def get_modality_suffix (modality):
    if modality[:3] == 'seg':
        return modality[3:]
    elif modality[:4] == 'hipp':
        return modality[4:]
    else:
        return modality



def get_final_spec():
    surf_spaces = []
    if 'cropT1w' in config['output_spaces']:
        surf_spaces.append('T1w')
    if 'corobl' in config['output_spaces']:
        surf_spaces.append('corobl')

    if len(config['hemi']) == 2:
        specs = expand(
            bids(root=root,datatype='surf',den='{density}',space='{space}',label='{autotop}', suffix='surfaces.spec', **config['subj_wildcards']),
                density=config['output_density'],
                space=surf_spaces,
                autotop=config['autotop_labels'],
                allow_missing=True)
    else:
         specs = expand(
            bids(root=root,datatype='surf',den='{density}',space='{space}',hemi='{hemi}',label='{autotop}',suffix='surfaces.spec', **config['subj_wildcards']),
                density=config['output_density'],
                space=surf_spaces,
                hemi=config['hemi'],
                autotop=config['autotop_labels'],
                allow_missing=True)
    return specs

def get_final_subfields():
    return expand(
        bids(
                root=root,datatype='seg',
                desc='subfields',suffix='dseg.nii.gz',
                space='{space}',hemi='{hemi}', **config['subj_wildcards']),
            hemi=config['hemi'],
            space=config['output_spaces'],
            allow_missing=True)

def get_final_coords():
    if 'laplace' in config['laminar_coords_method']:
        desc_io = 'laplace'
    elif 'equivolume' in config['laminar_coords_method']:
        desc_io = 'equivol'

    coords = []
    #compute all laplace coords by default (incl IO)
    coords.extend(
                expand(
                    bids(
                        root=root,datatype='seg',dir='{dir}',suffix='coords.nii.gz', desc='{desc}',space='{space}',hemi='{hemi}', label='{autotop}', **config['subj_wildcards']),
                            desc='laplace',
                            dir=['AP','PD','IO'],
                            autotop=config['autotop_labels'],
                            hemi=config['hemi'],
                            space=config['output_spaces'],
                            allow_missing=True))
    coords.extend(
                expand(
                    bids(
                        root=root,datatype='seg',dir='{dir}',suffix='coords.nii.gz', desc='{desc}',space='{space}',hemi='{hemi}', label='hipp', **config['subj_wildcards']),
                            desc=[desc_io],
                            dir=['IO'],
                            hemi=config['hemi'],
                            space=config['output_spaces'],
                            allow_missing=True))
    return coords


def get_final_transforms():
    if 'cropT1w' in  config['output_spaces']: 
        output_ref = 'T1w'
    else:
        output_ref = 'corobl'

    xfms = []

    xfms.extend(expand(
        bids(
                root=root,
                datatype='seg',
                **config['subj_wildcards'],
                label='{autotop}',
                suffix='xfm.nii.gz',
                hemi='{hemi}',
                from_='{space}',
                to='unfold',
                mode='image'),
            space=output_ref,
            autotop=config['autotop_labels'],
            hemi=config['hemi'],
            allow_missing=True))

    xfms.extend(expand(
        bids(
                root=root,
                datatype='seg',
                **config['subj_wildcards'],
                label='{autotop}',
                suffix='xfm.nii.gz',
                hemi='{hemi}',
                from_='unfold',
                to='{space}',
                mode='image'),
            space=output_ref,
            autotop=config['autotop_labels'],
            hemi=config['hemi'],
            allow_missing=True))
    xfms.extend(expand(
        bids(
                root=root,
                datatype='seg',
                **config['subj_wildcards'],
                label='{autotop}',
                suffix='xfm.nii.gz',
                hemi='{hemi}',
                from_='{space}',
                to='unfold',
                mode='image'),
            space=output_ref,
            autotop=config['autotop_labels'],
            hemi=config['hemi'],
            allow_missing=True))
    xfms.extend(expand(
        bids(
                root=root,
                datatype='seg',
                **config['subj_wildcards'],
                label='{autotop}',
                suffix='refvol.nii.gz',
                space='unfold'),
            autotop=config['autotop_labels'],
            allow_missing=True))

    return xfms


def get_final_anat():
    anat = []
    
    if 'cropT1w' in config['output_spaces']:
        anat.extend(
            expand(
                bids(
                        root=root,
                        datatype='seg',
                        desc='preproc',
                        suffix='{modality_suffix}.nii.gz'.format(
                                modality_suffix=get_modality_suffix(config['modality'])),
                        space='{space}',
                        hemi='{hemi}',
                        **config['subj_wildcards']),
                    space=config['output_spaces'],
                    hemi=config['hemi'],
                    allow_missing=True))
    return anat


def get_final_qc():
    qc = []
    #right now can only do qc from cropT1w space 
    if 'cropT1w' in config['output_spaces']:

        qc.extend(
            expand(
                bids(
                        root=root,
                        datatype='qc',
                        suffix='regqc.png',
                        from_='subject', 
                        to=config['template'],
                        **config['subj_wildcards']),
                    allow_missing=True)
            )

        qc.extend(
            expand(
                bids(
                        root=root,
                        datatype='qc',
                        suffix='dseg.png',
                        desc='subfields',
                        space='cropT1w',
                        hemi='{hemi}',
                        **config['subj_wildcards']),
                    hemi=config['hemi'],
                    allow_missing=True)
            )
        qc.extend(
            expand(
                bids(
                        root=root,
                        datatype='qc',
                        suffix='midthickness.surf.png', 
                        den='{density}',
                        desc='subfields',
                        space='cropT1w',
                        hemi='{hemi}',
                        label='{autotop}',
                        **config['subj_wildcards']),
                    hemi=config['hemi'],
                    autotop=config['autotop_labels'],
                    density=config['output_density'],
                    allow_missing=True)
            ) 
        if len(config['hemi']) == 2:
            qc.extend(
                expand(
                    bids(
                            root=root,
                            datatype='qc',
                            desc='subfields',
                            suffix='volumes.png',
                            **config['subj_wildcards']),
                        allow_missing=True)
                )

        if (config['modality'] == 'T1w' ) or (config['modality'] == 'T2w'):
            qc.extend(
                expand(
                    bids(
                            root=root,
                            datatype='qc',
                            desc='unetf3d',
                            suffix='dice.tsv',
                            hemi='{hemi}',
                            **config['subj_wildcards']),
                        hemi=config['hemi'],
                        allow_missing=True)
                )
    return qc


def get_final_subj_output():
    subj_output = []
    subj_output.extend(get_final_spec())
    subj_output.extend(get_final_subfields())
    subj_output.extend(get_final_coords())
    subj_output.extend(get_final_transforms())
    subj_output.extend(get_final_anat())
    subj_output.extend(get_final_qc())
    return subj_output

   

def get_final_output():

    if config['keep_work']:
        subj_output = get_final_subj_output()
    else:
        subj_output = get_final_work_tar()

    final_output = []

    modality_suffix = get_modality_suffix(config['modality'])
    modality_key = get_modality_key(config['modality'])

    final_output.extend(
        expand(subj_output,
                modality_suffix=modality_suffix,
                subject=config['input_lists'][modality_key]['subject'],
                session=config['sessions'])
        )

    return final_output


rule copy_to_results:
    """ Generic rule for copying data from work to results"""
    input:
        os.path.join(work, "{file}"),
    output:
        os.path.join(root, "{file}"),  
    group: 'subj'
    shell: 'cp {input} {output}'


def get_final_work_tar():
    return bids(root=work,suffix='work.tar.gz',
                include_subject_dir=False,
                include_session_dir=False,
                **config['subj_wildcards'])


def get_work_dir(wildcards):
    folder_with_file = expand(bids(root=work,**config['subj_wildcards']),**wildcards)
    folder_without_file = os.path.dirname(folder_with_file[0])
    return folder_without_file


rule archive_work_after_final:
    input: get_final_subj_output()
    params:
        work_dir = get_work_dir
    output: get_final_work_tar()
    group: 'subj'
    shell: 'tar -cvzf {output} {params.work_dir}; '
           'if [ $? -le 1 ]; then ' #exit code 0 or 1 is acceptable (2 is fatal)
           '  rm -rf {params.work_dir}; '
           'else exit 1; '
           'fi'


def get_input_for_shape_inject(wildcards):
    if config['modality'] == 'cropseg':
        seg = bids(root=work,datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='cropped',space='corobl',hemi='{hemi}').format(**wildcards)
    elif get_modality_key(config['modality']) == 'seg':
        modality_suffix = get_modality_suffix(config['modality'])
        seg = bids(root=work,datatype='anat',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='cropped',space='corobl',hemi='{hemi}',from_='{modality_suffix}').format(
                    **wildcards, modality_suffix=modality_suffix),
    else:
        seg = bids(root=work,datatype='seg',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='nnunet',space='corobl',hemi='{hemi}').format(**wildcards)
    return seg

def get_labels_for_laplace(wildcards):
    if config['skip_inject_template_labels']:
        seg = get_input_for_shape_inject(wildcards)
    else:
        seg = bids(root=work,datatype='seg',**config['subj_wildcards'],suffix='dseg.nii.gz',desc='postproc',space='corobl',hemi='{hemi}').format(**wildcards)
    return seg

