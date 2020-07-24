#!/usr/bin/env python3
import argparse
import os
import subprocess
import nibabel
import numpy
from bids import BIDSLayout
from bids.layout import Config
from glob import glob

__version__ = open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'version')).read()

def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d"%process.returncode)

parser = argparse.ArgumentParser(description='Hippocampal AutoTop BIDS App')
parser.add_argument('bids_dir', help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('output_dir', help='The directory where the output files '
                    'should be stored. If you are running group level analysis '
                    'this folder should be prepopulated with the results of the'
                    'participant level analysis.')
parser.add_argument('analysis_level', help='Level of the analysis that will be performed. '
                    'Multiple participant level analyses can be run independently '
                    '(in parallel) using the same output_dir.',
                    choices=['participant', 'group'])
parser.add_argument('--participant_label', help='The label(s) of the participant(s) that should be analyzed. The label '
                   'corresponds to sub-<participant_label> from the BIDS spec '
                   '(so it does not include "sub-"). If this parameter is not '
                   'provided all subjects should be analyzed. Multiple '
                   'participants can be specified with a space separated list.',
                   nargs="+")
parser.add_argument('--suffix', help='Only use images with the specified suffix entity in the filename',
                    default='T2w')
parser.add_argument('--acq', help='Only use images with the specified acq entity in the filename')
parser.add_argument('--run', help='Only use images with the specified run entity in the filename')
parser.add_argument('--search', help='Wildcard search term to locate in image, use when multiple T2w images match for a subject')

#parser.add_argument('--skip_bids_validator', help='Whether or not to perform BIDS dataset validation',
#                   action='store_true')
parser.add_argument('-v', '--version', action='version',
                    version='Hippocampal AutoTop BIDS App version {}'.format(__version__))


args = parser.parse_args()

#if not args.skip_bids_validator:
#    run('bids-validator %s'%args.bids_dir)

#default search term is by suffix
search_terms = {'suffix': args.suffix}

#add optional search terms
if args.acq is not None:
    search_terms['acquisition'] = args.acq
if args.run is not None:
    search_terms['run'] = args.run



#start by getting bids layout from pybids
layout = BIDSLayout(args.bids_dir,derivatives=False,validate=False,index_metadata=False)

subjects_to_analyze = []
# only for a subset of subjects
if args.participant_label:
    subjects_to_analyze = args.participant_label
# for all subjects
else:
    subject_dirs = glob(os.path.join(args.bids_dir, "sub-*"))
    subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]

#print(subjects_to_analyze)
#subjects_to_analyze = subjects_to_analyze.sort()
#print(subjects_to_analyze)

in_imgs=[]
out_dirs=[]

# running participant level
if args.analysis_level == "participant":

    for subject in sorted(subjects_to_analyze):

        t2w_list = layout.get(subject=subject,extension='nii.gz',return_type='filename',**search_terms)
        if args.search and len(t2w_list)>0:
            t2w_list = [t2w for t2w in t2w_list if args.search in t2w ]

        if len(t2w_list) == 0:
            print(f'WARNING: Unable to find any matching images for sub-{subject}!')
            continue
        elif len(t2w_list) >1:
            print(f'WARNING: Found multiple matching T2w files for sub-{subject}, please use search filters to narrow down to one')
            print(t2w_list)
            continue
        elif len(t2w_list) == 1:
#            print(f'Found image file for sub-{subject}')
#            print(t2w_list)
            in_imgs = in_imgs + t2w_list #in_img
            out_dirs = out_dirs + [os.path.join(args.output_dir,f'sub-{subject}')]
            
            
    for in_img, out_dir in zip(in_imgs,out_dirs):
        print(f'/src/mcr_v97/run_singleSubject.sh /opt/mcr/v97 {in_img} {out_dir}')

# running group level
elif args.analysis_level == "group":
    print('insert myousif report generation here!')

