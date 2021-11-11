import numpy as np
import nibabel as nib
from astropy.convolution import convolve as nan_convolve
from scipy.interpolate import NearestNDInterpolator
import skfmm


logfile = open(snakemake.log[0], 'w')

# This script will get approximate Laplace_PD within the DG by fast marching edge-to-edge within small AP slices. Most distant voxels are used as endpoints. Start in the middle of the hippocmapus and then NN copy to next slice and rerun fast marching both directions (to ensure the march goes in the same direction between slices!). Note that this requires endpoints to be most distant from eachother in each slice, as in a 'U' shape, but this isn't always the case if there is more of a 'V' shape!

# test inputs
#lbl_nib = nib.load('test_T1w_dentate-noCA4/work/work/sub-01/seg_T1w/sub-01_hemi-R_space-corobl_desc-postproc_dseg.nii.gz')
#lbl = lbl_nib.get_fdata()
#gmlbl = [8]
#AP_nib = nib.load('test_T1w_dentate-noCA4/work/work/sub-01/seg_T1w/sub-01_dir-AP_hemi-R_space-corobl_desc-laplace_coords.nii.gz')
#AP = AP_nib.get_fdata()

# real inputs
lbl_nib = nib.load(snakemake.input.lbl)
lbl = lbl_nib.get_fdata()
gmlbl = snakemake.params.gm_labels
AP_nib = nib.load(snakemake.input.APcoords)
AP = AP_nib.get_fdata()
print(f'data loaded', file=logfile, flush=True)

# params
nslices = snakemake.params.nslices # 30 seems to work well for standard iamges
smooth_iters = snakemake.params.smooth_iters # still testing. 5 seems fine.
ap = np.linspace(0,1,nslices)
speed = np.zeros_like(lbl) + 0.01 # allow slow travel outside DG
speed[lbl==gmlbl] = 1
PD = np.zeros_like(AP)
PD[:] = np.nan

# start with middle AP slice
i = nslices//2
sl = np.logical_and(AP>ap[i-1], AP<ap[i])
APslice_mid = np.logical_and(sl,lbl==gmlbl)
[x,y,z] = np.where(APslice_mid)
med = np.median(x) # sagittally the middle
# march outward (both laterally and medially)
v = np.where(x==med)[0]
phi = np.ones_like(lbl)
phi[x[v], y[v], z[v]] = 0
forward = skfmm.travel_time(phi,speed)
out = np.zeros_like(forward)
out[APslice_mid] = forward[APslice_mid]
v = np.where(x<med)[0] # block out lateral side so we always start medial
out[x[v], y[v], z[v]] = 0
# now march lateral
[x,y,z] = np.where(out==np.max(out))
phi = np.ones_like(lbl)
phi[x[0],y[0],z[0]] = 0 # only the first point
forward = skfmm.travel_time(phi,speed)
PD[APslice_mid] = forward[APslice_mid]/np.max(forward[APslice_mid])
print(f'fastmarch for middle AP slice done', file=logfile, flush=True)

# move to new slices (towards anterior)
APslice_prev = APslice_mid
for ii in range(i-1,0,-1):
    sl = np.logical_and(AP>ap[ii-1], AP<ap[ii])
    APslice = np.logical_and(sl,lblgmlbl)
    if np.sum(APslice) > 0:
        # interpolate this slice
        [x,y,z] = np.where(APslice_prev)
        interp = NearestNDInterpolator(np.c_[x,y,z],PD[APslice_prev])
        [x,y,z] = np.where(APslice)
        out = interp(x,y,z)
        # march backward
        v = np.where(out==(np.max(out)))[0] # most lateral
        phi = np.ones_like(lbl)
        phi[x[v], y[v], z[v]] = 0
        forward = skfmm.travel_time(phi,speed)
        out = forward[APslice]
        # march forward
        v = np.where(out==(np.max(out)))[0] # most lateral
        phi = np.ones_like(lbl)
        phi[x[v], y[v], z[v]] = 0
        forward = skfmm.travel_time(phi,speed)
        PD[APslice] = forward[APslice]/np.max(forward[APslice])
        print('fastmarch for AP slice {ii} done', file=logfile, flush=True)
    else:
        print('skipping AP slice {ii} not enough voxels', file=logfile, flush=True)
    APslice_prev = APslice

# move to new slices (towards posterior)
APslice_prev = APslice_mid
for ii in range(i+1,nslices,1):
    sl = np.logical_and(AP>ap[ii-1], AP<ap[ii])
    APslice = np.logical_and(sl,lblgmlbl)
    if np.sum(APslice) > 0:
        # interpolate this slice
        [x,y,z] = np.where(APslice_prev)
        interp = NearestNDInterpolator(np.c_[x,y,z],PD[APslice_prev])
        [x,y,z] = np.where(APslice)
        out = interp(x,y,z)
        # march backward
        v = np.where(out==(np.max(out)))[0] # most lateral
        phi = np.ones_like(lbl)
        phi[x[v], y[v], z[v]] = 0
        forward = skfmm.travel_time(phi,speed)
        out = forward[APslice]
        # march forward
        v = np.where(out==(np.max(out)))[0] # most lateral
        phi = np.ones_like(lbl)
        phi[x[v], y[v], z[v]] = 0
        forward = skfmm.travel_time(phi,speed)
        PD[APslice] = forward[APslice]/np.max(forward[APslice])
        print(f'fastmarch for AP slice {ii} done', file=logfile, flush=True)
    else:
        print(f'skipping AP slice {ii} not enough voxels', file=logfile, flush=True)
    APslice_prev = APslice

# smooth to clean up space between slices
hl=np.zeros([3,3,3])
hl[1,:,:] = 1
hl[:,1,:] = 1
hl[:,:,1] = 1
hl[1,1,1] = 0
hl = hl/np.sum(hl)
PD_smooth = PD
for n in range(smooth_iters):
    PD_smooth = nan_convolve(PD_smooth,hl,preserve_nan=True)

print(f'smoothing done', file=logfile, flush=True)

sv = nib.Nifti1Image(PD_smooth,AP_nib.affine,AP_nib.header)
nib.save(sv,snakemake.output.coords_pd)
logfile.close()


