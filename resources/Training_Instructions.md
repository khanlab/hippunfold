This document contains step-by-step instructions for retraining or finetuning UNet for hippocampal segmentation. 

## 1) Run existing data through a previous model
There's a chance that one of the existing models will perform well on your new dataset, if the data is similar enough. If performance is very good then no further fine-tuning is needed. If performance is poor on some samples, they can be manually corrected or else fully manually segmented. In either case, running the full pipeline end-to-end should produce images that are CorObl, which is the space that segmentations for training should be in.

see `hippunfold -h`

## 2) Collect training images and segmentations
All training data should be manually inspected, and once the quality is good the CorObl image (eg. `outputdirectory/subjectID/hemi-L/img.nii`) and corresponding segmentation image (eg. `outputdirectory/subjectID/hemi-L/niftynet_lbl.nii` or a manually generated segmentation image) can be copied into a new clean directory (eg. `mynewdataset/training/`). Each new subject in the training directory should have a unique subjectID as a prefix, and either `_img` or `_lbl` for images and segmentations, repsectively.

For example `ls mynewdataset/training` should produce something like this:
```
sub-001_img.nii.gz
sub-001_lbl.nii.gz
sub-002_img.nii.gz
sub-002_lbl.nii.gz
```
It is also possible to fine-tune on only a subset of subjects (for example, only those that produced good performance on the first pass).

## 3) Fine-tune an existing model, or train one from scratch
Once you have populated your training data directory, you may train your model or fine-tune an existing model using `fineTune_UNet.sh`. This is a compute-intensive process. This can be run on a CPU, but it is recommended that you run on GPU with sufficient GPU memory (current models were trained on 8xV100 GPU nodes). By default, `fineTune_UNet.sh` will run 100k iterations which should take <\24h with these parameters. 

For example:
`singularity exec --nv hippocampal_autotop_latest.sif bash /src/resources/fineTune_UNet.sh mynewdataset/training mynewdataset/newCNNmodel`
(omit `--nv` if no GPU is available)

This will perform 100k training iterations with data augmentation and using the same parameters as previous work. Training and validation progress can be viewed using tensorboard (eg. `tensorboard --logdir mynewdataset/newCNNmodel/models`). Once training is complete, inference will be performed on the remaining test data, which can then be inspected for quality. Further training iterations can be run using the same command as above (specifying the same output directory), or a new model can be trained using the same data by specifying a different directory. 

If you know what you are doing, you can open `mynewdataset/newCNNmodel/config.ini` and modify parameters before running additional training. 

## 4) Incremental learning
If your dataset is very large, you may fine-tune on only a subset of new samples. In that case, you can re-run steps 1-3 which should now produce more good quality segmentations for use in further training. 

## 5) Share trained models and/or data
Please consider sharing your data and/or trained models to improve generalizability to future studies.

