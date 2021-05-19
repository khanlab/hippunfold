import nibabel as nib
import numpy as np
import scipy.ndimage

# this function solves the Laplace equation for Anterior-Posterior, Proximal-distal, and Inner-Outer axes of the hippocamps
maxiters = [10000,5000,1000]
nii_labelmap = snakemake.input.lbl

### testing only ###
# nii_labelmap = 'test/hemi-L/labelmap-postProcess.nii.gz'
nii_ap = snakemake.output.coords_ap
nii_pd = snakemake.output.coords_pd
nii_io = snakemake.output.coords_io

# TODO: initialize with fast marching for faster solver
# TODO: make iters a parameter or choose stopping contitions or both (these values should be fine for 0.3mm resolution)

# load labelmap
labelmap = nib.load(nii_labelmap)
lbl = labelmap.get_fdata()
idxgm = np.zeros(lbl.shape)
idxgm[lbl==1] = 1
idxgm[lbl==8] = 1

# AP solution
source = np.zeros(lbl.shape)
source[lbl==5] = 1
sink = np.zeros(lbl.shape)
sink[lbl==6] = 1
AP = laplace_iters(idxgm,source,sink,iters=maxiters[0])
coordsAP = nib.Nifti1Image(AP,labelmap.affine,labelmap.header)
nib.save(coordsAP,nii_ap)

# PD solution
source = np.zeros(lbl.shape)
source[lbl==3] = 1
sink = np.zeros(lbl.shape)
sink[lbl==8] = 1
PD = laplace_iters(idxgm,source,sink,iters=maxiters[1])
coordsPD = nib.Nifti1Image(PD,labelmap.affine,labelmap.header)
nib.save(coordsPD,nii_pd)

# IO solution
source = np.zeros(lbl.shape)
source[lbl==2] = 1
source[lbl==4] = 1
source[lbl==7] = 1
sink = np.zeros(lbl.shape)
sink[lbl==0] = 1
IO = laplace_iters(idxgm,source,sink,iters=maxiters[2])
coordsIO = nib.Nifti1Image(IO,labelmap.affine,labelmap.header)
nib.save(coordsIO,nii_io)

def coords = laplace_iters(idxgm,source,sink,iters=100):
	# set up filter (18NN)
	hl=np.zeros([3,3,3])
	hl[1,:,:] = 1
	hl[:,1,:] = 1
	hl[:,:,1] = 1
	hl[1,1,1] =
	hl = hl/np.sum(hl)
	bg = 1-idxgm
	bg[source==1] = 0
	bg[sink==1] = 0
	bgCorrection = scipy.ndimage.convolve(bg,hl)
	coords = np.zeros(lbl.shape)
	coords[idxgm==1] = 0.5
	coords[sink==1] = 1
	for i = range(100):
		coords = scipy.ndimage.convolve(AP,hl)
		coords = AP/bgCorrection
		coords[bg==1] = 0
		coords[source==1] = 0
		coords[sink==1] = 1
	coords[idxgm==0] = 0
return coords
