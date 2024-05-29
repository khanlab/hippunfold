# Template-base segmentation

Template-based segmentation can be used with the `--use-template-seg` flag instead of a deep neural network for tissue class segmentation prior to unfolding. This is the recommended workflow for non-human data.

**Advantages:**
- Relatively robust to image quality
- Can be used without UNet training (which requires many manually segmented samples)
- Precision can be adjusted with the `--inject_template_smoothing_factor` and `--rigid-reg-template` flags

**Disadvantages:**
- Doesn't account well for interindividual differences in folding patterns (not an issue in most non-human species where the hippocampus is relatively smooth)
- Can still fail due to registration errors

This is meant as an alternative to UNet-based tissue segmentation when only one or a few manually segmented training samples are available, which are not sufficient to train a UNet model. However, failures can still occur during this registration due to differences in image contrast and/or quality compared to the template. Adjusting the optional parameters is sufficient to solve this in most cases, but if not then you can manually register or segment your hippocampal images and then run with the `--modality segT1w`, `--modality segT2w`, or `--modality cropseg` (the latter doesn't perform any additional cropping its its recommended you crop your segmentations to improve HippUnfold processing time). 


If you have a unique template hippocmapal segmentation (e.g. from another species or special population) then please consider making it available to other HippUnfold users! We will be happy to include it, and any associated references, if you raise a [git issue](https://github.com/khanlab/hippunfold/issues). 
