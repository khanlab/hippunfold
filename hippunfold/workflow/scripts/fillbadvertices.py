import nibabel as nib
import numpy as np
from scipy.stats import zscore
import copy

SDthreshold = snakemake.params.threshold
iters = snakemake.params.dist

gii = nib.load(snakemake.input.gii)
varr = gii.get_arrays_from_intent("NIFTI_INTENT_POINTSET")[0]
V = varr.data
farr = gii.get_arrays_from_intent("NIFTI_INTENT_TRIANGLE")[0]
F = farr.data


# find local outliers by smoothing and then substracting from original
# https://github.com/MICA-MNI/hippomaps/blob/master/hippomaps/utils.py
def avg_neighbours(F, cdat, n):
    frows = np.where(F == n)[0]
    v = np.unique(F[frows, :])
    cdat = np.reshape(cdat, (len(cdat), -1))
    out = np.nanmean(cdat[v, :], 0)
    return out


def surfdat_smooth(F, cdata, iters=1):
    sz = cdata.shape
    cdata = cdata.reshape(cdata.shape[0], -1)
    cdata_smooth = copy.deepcopy(cdata)
    for i in range(iters):
        for n in range(len(cdata)):
            cdata_smooth[n, :] = avg_neighbours(F, cdata, n)
        cdata = copy.deepcopy(cdata_smooth)
    return cdata_smooth.reshape(sz)


Vsmooth = surfdat_smooth(F, V, iters=iters)
Vdiffs = V - Vsmooth
Vdists = np.sqrt((Vdiffs[:, 0]) ** 2 + (Vdiffs[:, 1]) ** 2 + (Vdiffs[:, 2]) ** 2)
Vzscored = zscore(Vdists)
outliers = (Vzscored > SDthreshold) | (Vzscored < -SDthreshold)
V[outliers, :] = np.nan


# most nans should be just isolated points, but in case there is an island of nans this will erode it, replacing with decent (but not perfect) guesses of where vertices should be
while np.isnan(np.sum(V)):
    # index of vertices containing nan
    i = np.where(np.isnan(V))
    ii = np.unique(i[0])
    # replace with the nanmean of neighbouring vertices
    newV = V
    for n in ii:
        f = np.where(F == n)
        v = F[f[0]]
        vv = np.unique(v)
        newV[n, :] = np.nanmean(V[vv, :], 0)
    V = newV

nib.save(gii, snakemake.output.gii)
