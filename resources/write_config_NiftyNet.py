# write config file for NiftyNet highres3dnet network training using these default parameters. First argument is the training directory, second argument is the output directory, third argument specifies number of iterations, and fourht (optional) an existing model to bootstrap.

import sys
import configparser
import os
import shutil
import glob

# input arguments

trainingdir = sys.argv[1] #'../training_data_b1000/'
newmodeldir = sys.argv[2] #'testmodel'
iterations = sys.argv[3] #'testmodel'
try:
	os.mkdir(newmodeldir)
except:
	print('output directory already exists')

# copy over bootstrapped CNN model
if len(sys.argv) == 5:
	bootstrapmodel = sys.argv[4]
	shutil.copytree(bootstrapmodel + '/models/', newmodeldir + '/models/')
	start_iter = '-1'
	# add existing iterations to max iterations
	fn = glob.glob(newmodeldir + '/models/*.index')[0]
	i = fn.find('ckpt-')[0]
	fn = fn[i+5:-6]
	iterations = int(iterations) + int(fn)
elif len(sys.argv) == 4:
	bootstrapmodel = '' # Optional. To resume training a timed out model, specify it here.
	start_iter = '0'
else:
	print('Error wrong number of input arguments')

# write config file with default values
config = configparser.ConfigParser()

config['IMG'] = {'path_to_search': trainingdir,
		'filename_contains': 'img',
             	'spatial_window_size': '(64, 64, 64)',
		'interp_order': '1',
		'pixdim': '(0.3, 0.3, 0.3)',
		'axcodes': '(R, A, S)'}
config['LBL'] = {'path_to_search': trainingdir,
                'filename_contains': 'lbl',
                'spatial_window_size': '(64, 64, 64)',
                'interp_order': '1',
                'pixdim': '(0.3, 0.3, 0.3)',
                'axcodes': '(R, A, S)'}

config['SYSTEM'] = {'cuda_devices': '""',
		'model_dir': newmodeldir}

config['NETWORK'] = {'name': 'highres3dnet_large',
		'batch_size': '1',
		'activation_function': 'relu',
		'volume_padding_size': '0',
		'normalisation': 'True',
		'foreground_type': 'mean_plus',
		'cutoff': '(0.001, 0.999)'}

config['TRAINING'] = {'sample_per_volume': '5',
		'lr': '0.001',
		'loss_type': 'Dice',
		'starting_iter': start_iter,
		'save_every_n': '1000',
		'tensorboard_every_n': '100',
		'max_iter': iterations,
		'validation_every_n': '100',
		'exclude_fraction_for_validation': '0.2',
		'exclude_fraction_for_inference': '0.2',
		'rotation_angle': '(-10.0,10.0)',
		'random_flipping_axes': '0',
		'do_elastic_deformation': 'True',
		'num_ctrl_points': '4',
		'deformation_sigma': '15',
		'proportion_to_deform': '0.75',
		'bias_field_range': '(-5.0,5.0)',
		'bf_order': '3'}


config['INFERENCE'] = {'border': '(16,16,16)',
		'inference_iter': '-1',
		'save_seg_dir': newmodeldir + '/parcellation_output',
		'output_interp_order': '0'}

config['SEGMENTATION'] = {'image': 'IMG',
		'label': 'LBL',
		'label_normalisation': 'False',
		'output_prob': 'False',
		'num_classes': '9'}

config['EVALUATION'] = {'save_csv_dir': newmodeldir + '/eval',
		'evaluations': 'dice,average_distance'}

with open(newmodeldir + '/config.ini', 'w') as configfile:
  config.write(configfile)

