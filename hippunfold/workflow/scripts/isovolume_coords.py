import os
import nighres
import nibabel
from shutil import copyfile

logfile = open(snakemake.log[0], 'w')


## first load labelmap and binarize one edge that includes the inner surf and one edge the outer
lbl_nib = nib.load(snakemake.input.lbl)
lbl = lbl_nib.get_fdata()
print(f'labelmap loaded', file=logfile, flush=True)

outfile = snakemake.output.coords
path,filename = os.path.split(outfile)

source = np.zeros(lbl.shape)
for i in snakemake.params.src_labels:
    source[lbl==i] = 1
bin_nib = nib.Nifti1Image(source,lbl_nib.affine,lbl_nib.header)
nib.save(bin_nib,path + '/innersurf_bin.nii.gz')
print(f'inner edge binarized', file=logfile, flush=True)

sink = np.zeros(lbl.shape)
sink[lbl>0] = 1
bin_nib = nib.Nifti1Image(sink,lbl_nib.affine,lbl_nib.header)
nib.save(bin_nib,path + '/outsurf_bin.nii.gz')
print(f'outer edge binarized', file=logfile, flush=True)


## convert binarized edges to levelset surfaces
nighres.surface.probability_to_levelset(path + '/innersurf_bin.nii.gz', save_data=True, output_dir=path)
nighres.surface.probability_to_levelset(path + '/outersurf_bin.nii.gz', save_data=True, output_dir=path)
print(f'binarized files to levelsets complete', file=logfile, flush=True)

## get isovolume solution
nighres.laminar.volumetric_layering(path + '/innersurf_bin_p2l-surf.nii.gz', 'outersurf_bin_p2l-surf.nii.gz', save_data=True, output_dir=path)
print(f'levelset to ', file=logfile, flush=True)

## copy to output (volumetric_layering produces multiple files so we need to specify the right one)
copyfile(path + 'innersurf_bin_p2l-surf_layering-depth.nii.gz', snakemake.output.coords)
