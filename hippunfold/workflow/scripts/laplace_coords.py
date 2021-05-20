import nibabel as nib
import numpy as np
import scipy.ndimage


def laplace_iters(idxgm,source,sink,iters=100):
    # set up filter (18NN)
    hl=np.zeros([3,3,3])
    hl[1,:,:] = 1
    hl[:,1,:] = 1
    hl[:,:,1] = 1
    hl[1,1,1] = 0
    hl = hl/np.sum(hl)
    bg = 1-idxgm
    bg[source==1] = 0
    bg[sink==1] = 0
    bgCorrection = scipy.ndimage.convolve(bg,hl)
    coords = np.zeros(lbl.shape)
    coords[idxgm==1] = 0.5
    coords[sink==1] = 1
        
    for i in range(iters): 
        coords = scipy.ndimage.convolve(coords,hl)
#        coords = coords/bgCorrection  #this is causing a divide by zero
        coords[bg==1] = 0
        coords[source==1] = 0
        coords[sink==1] = 1

    coords[idxgm==0] = 0

    return coords


# this function solves the Laplace equation for Anterior-Posterior, Proximal-distal, and Inner-Outer axes of the hippocamps
nii_labelmap = snakemake.input.lbl

iters_ap = snakemake.params.iters_ap
iters_pd = snakemake.params.iters_pd
iters_io = snakemake.params.iters_io

nii_ap = snakemake.output.coords_ap
nii_pd = snakemake.output.coords_pd
nii_io = snakemake.output.coords_io

# TODO: initialize with fast marching for faster solver
# TODO: make iters a parameter or choose stopping contitions or both (these values should be fine for 0.3mm resolution)

# load labelmap
print('loading labelmap')
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
print('performing laplace on AP')
AP = laplace_iters(idxgm,source,sink,iters=iters_ap)
print('done performing laplace on AP')
coordsAP = nib.Nifti1Image(AP,labelmap.affine,labelmap.header)
nib.save(coordsAP,nii_ap)

# PD solution
source = np.zeros(lbl.shape)
source[lbl==3] = 1
sink = np.zeros(lbl.shape)
sink[lbl==8] = 1
PD = laplace_iters(idxgm,source,sink,iters=iters_pd)
coordsPD = nib.Nifti1Image(PD,labelmap.affine,labelmap.header)
nib.save(coordsPD,nii_pd)

# IO solution
source = np.zeros(lbl.shape)
source[lbl==2] = 1
source[lbl==4] = 1
source[lbl==7] = 1
sink = np.zeros(lbl.shape)
sink[lbl==0] = 1
IO = laplace_iters(idxgm,source,sink,iters=iters_io)
coordsIO = nib.Nifti1Image(IO,labelmap.affine,labelmap.header)
nib.save(coordsIO,nii_io)


