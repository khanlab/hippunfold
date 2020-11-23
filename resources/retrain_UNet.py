import configparser
import os
import shutil



# input arguments
# TODO: make this a proper function
trainingdir = '../training_data_b1000/'
newmodeldir = 'testmodel'
bootstrapmodel = '' # Optional. To resume training a timed out model, specify it here.



os.mkdir(newmodeldir)
if bootstrapmodel != '':
	shutil.copytree(bootstrapmodel + '/models/', newmodeldir + '/models/')
	start_iter = '-1'
else:
	start_iter = '0'





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
		'max_iter': '500000',
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
		'label': 'lbl',
		'label_normalisation': 'False',
		'output_prob': 'False',
		'num_classes': '9'}

config['EVALUATION'] = {'save_csv_dir': newmodeldir + '/eval',
		'evaluations': 'dice,average_distance'}

with open(newmodeldir + '/config.ini', 'w') as configfile:
  config.write(configfile)


# TO BE RUN IN BASH
# requires niftynet
net_segment -c newmodeldir/config.ini train
net_segment -c newmodeldir/config.ini inference
net_segment -c newmodeldir/config.ini evaluation


