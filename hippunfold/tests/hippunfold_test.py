import os
import tempfile


def test_dry_runs(script_runner):
    os.environ['HIPPUNFOLD_CACHE_DIR'] = os.path.join(os.getcwd(),'test_data/fake_models')
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_singleT2w',output_dir,'participant','-np')
        assert ret.success

    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_multiT2w',output_dir,'participant','-np')
        assert ret.success

    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_T1w',output_dir,'participant','-np','--modality','T1w')
        assert ret.success
        
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_hippb500',output_dir,'participant','-np','--modality','hippb500')
        assert ret.success
        
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_T1w_longitudinal',output_dir,'participant','-np','--modality','T1w')
        assert ret.success

    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_singleT2w_longitudinal',output_dir,'participant','-np')

    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_segT2w',output_dir,'participant','-np','--modality','segT2w')
        assert ret.success

    #test case for cropseg, uses --path_cropseg instead of bids, since pybids (or at least how we use it in snakebids) doesn't like hemi wildcards
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', '-',output_dir,'participant','-np','--modality','cropseg','--path_cropseg','test_data/data_cropseg/sub-{subject}_hemi-{hemi}_dseg.nii.gz')
        assert ret.success

