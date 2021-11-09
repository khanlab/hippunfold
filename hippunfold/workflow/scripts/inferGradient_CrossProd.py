import numpy as np
import nibabel as nib
from astropy.convolution import convolve as nan_convolve
from scipy.interpolate import NearestNDInterpolator

logfile = open(snakemake.log[0], 'w')

# test inputs
#lbl_nib = nib.load('test_T1w_dentate-noCA4/work/work/sub-01/seg_T1w/sub-01_hemi-R_space-corobl_desc-postproc_dseg.nii.gz')
#lbl = lbl_nib.get_fdata()
#gmlbl = [8]
#AP_nib = nib.load('test_T1w_dentate-noCA4/work/work/sub-01/seg_T1w/sub-01_dir-AP_hemi-R_space-corobl_desc-laplace_coords.nii.gz')
#AP = AP_nib.get_fdata()
#APsrc = [5]
#APsnk = [6]
#IO_nib = nib.load('test_T1w_dentate-noCA4/work/work/sub-01/seg_T1w/sub-01_dir-IO_hemi-R_space-corobl_desc-dentate_coords.nii.gz')
#IO = IO_nib.get_fdata()
#IOsrc = [1]
#IOsnk = [2, 4, 7]

# real inputs
lbl_nib = nib.load(snakemake.input.lbl)
lbl = lbl_nib.get_fdata()
gmlbl = snakemake.params.lblgm
AP_nib = nib.load(snakemake.input.APcoords)
AP = AP_nib.get_fdata()
APsrc = snakemake.params.APsrc_labels
APsnk = snakemake.params.APsink_labels
IO_nib = nib.load(snakemake.input.IOcoords)
IO = IO_nib.get_fdata()
IOsrc = snakemake.params.IOsrc_labels
IOsnk = snakemake.params.IOsink_labels

print(f'data loaded', file=logfile, flush=True)

# pad with src/snk and NN interpolation so gradients dont jump at boundary
idxgm = np.zeros(lbl.shape)
for i in gmlbl:
    idxgm[lbl==i] = 1

source = np.zeros(lbl.shape)
for i in APsrc:
    source[lbl==i] = 1

sink = np.zeros(lbl.shape)
for i in APsnk:
    sink[lbl==i] = 1

AP[idxgm==0] = np.nan
AP[source==1] = 0
AP[sink==1] = 1

[x,y,z] = np.where(np.invert(np.isnan(AP)))
interp = NearestNDInterpolator(np.c_[x,y,z],AP[np.invert(np.isnan(AP))])
[qx,qy,qz] = np.where(np.isnan(AP))
fill = interp(qx,qy,qz)
AP[np.isnan(AP)] = fill

source = np.zeros(lbl.shape)
for i in IOsrc:
    source[lbl==i] = 1

sink = np.zeros(lbl.shape)
for i in IOsnk:
    sink[lbl==i] = 1

IO[idxgm==0] = np.nan
IO[source==1] = 0
IO[sink==1] = 1

[x,y,z] = np.where(np.invert(np.isnan(IO)))
interp = NearestNDInterpolator(np.c_[x,y,z],IO[np.invert(np.isnan(IO))])
[qx,qy,qz] = np.where(np.isnan(IO))
fill = interp(qx,qy,qz)
IO[np.isnan(IO)] = fill

# now take the cross product of the gradients
[IOfx,IOfy,IOfz] = np.gradient(IO)
[APfx,APfy,APfz] = np.gradient(AP)
PDfx = np.zeros_like(IOfx)
PDfy = np.zeros_like(IOfx)
PDfz = np.zeros_like(IOfx)
[x,y,z] = np.where(idxgm)
for i in range(len(x)):
    vox = x[i], y[i], z[i]
    c = np.cross([IOfx[vox], IOfy[vox], IOfz[vox]], [APfx[vox], APfy[vox], APfz[vox]])
    c[np.isnan(c)] = 0
    PDfx[vox] = c[0]
    PDfy[vox] = c[1]
    PDfz[vox] = c[2]

# inverse (integral) gradient (code from https://www.mathworks.com/matlabcentral/fileexchange/9734-inverse-integrated-gradient)
# that was too hard to translate (sparse matrix maths not supported)
# instead: home-brew integration
dxintegral = np.cumsum(PDfx,0)
dyintegral = np.cumsum(PDfy,1)
dzintegral = np.cumsum(PDfz,2)
intPD = np.zeros(PDfx.shape)
[x,y,z] = np.where(idxgm)
for i in range(len(x)):
    vox = x[i], y[i], z[i]
    intPD[vox] = dxintegral[vox] + dyintegral[vox] + dzintegral[vox]

print(f'inference complete', file=logfile, flush=True)

# smooth
out_smooth = intPD
hl=np.ones([3,3,3])
hl = hl/np.sum(hl)
out_smooth[idxgm==0] = np.nan
out_smooth = nan_convolve(out_smooth,hl,preserve_nan=True)
out_smooth[idxgm==0] = 0.

# rescale (per AP slice)
nslices = 20
out_norm = np.zeros(intPD.shape)
ap = np.linspace(0,1,nslices)
for i in range(len(ap-1)):
    sl = np.logical_and(AP<ap[i-1], AP>ap[i])
    slice = np.logical_and(sl, idxgm)
    if np.sum(slice) > 10:
        s = out_smooth[slice]
        s = s-np.nanmin(s)
        s = s/np.nanmax(s)
        out_norm[slice] = s
    else:
        out_norm[slice] = out_smooth[slice]

# smooth again?
out_smooth = out_norm
hl=np.ones([3,3,3])
hl = hl/np.sum(hl)
out_smooth[idxgm==0] = np.nan
out_smooth = nan_convolve(out_smooth,hl,preserve_nan=True)
out_smooth[idxgm==0] = 0.

# save
out = out_smooth
sv = nib.Nifti1Image(out,AP_nib.affine,AP_nib.header)
nib.save(sv,snakemake.output.coords_pd)
#nib.save(sv,'test')

print(f'smoothed, scaled, saved', file=logfile, flush=True)
logfile.close()
