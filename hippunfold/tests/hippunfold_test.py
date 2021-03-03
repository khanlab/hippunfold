import os
import tempfile


os.environ['HIPPUNFOLD_CACHE_DIR'] = os.path.join(os.getcwd(),'test_data/fake_models')

def test_dry_runs_T2w(script_runner):
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_singleT2w',output_dir,'participant','-np')
        assert ret.success

def test_dry_run_multiT2w(script_runner):
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_multiT2w',output_dir,'participant','-np')
        assert ret.success

def test_dry_run_T1w(script_runner):
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_T1w',output_dir,'participant','-np','--modality','T1w')
        assert ret.success
        
def test_dry_run_hippb500(script_runner):
    with tempfile.TemporaryDirectory() as output_dir:
        ret = script_runner.run('hippunfold', 'test_data/bids_hippb500',output_dir,'participant','-np','--modality','hippb500')
        assert ret.success
                
        
        
