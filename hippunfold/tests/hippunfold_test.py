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
        ret = script_runner.run('hippunfold', 'test_data/bids_T1w_longitudinal',output_dir,'participant','-np','--modality','T1w')
        assert ret.success

    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_singleT2w_longitudinal',output_dir,'participant','-np')
        assert ret.success


