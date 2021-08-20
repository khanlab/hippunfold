import os
import nighres
from shutil import copyfile

logfile = open(snakemake.log[0], 'w')
print(f'start', file=logfile, flush=True)

outfile = snakemake.output.coords
path,filename = os.path.split(outfile)

## convert binarized edges to levelset surfaces
nighres.surface.probability_to_levelset(snakemake.input.innerbin, save_data=True, output_dir=path)
nighres.surface.probability_to_levelset(snakemake.input.outerbin, save_data=True, output_dir=path)
print(f'binarized files to levelsets complete', file=logfile, flush=True)

## get isovolume solution
inbin = snakemake.input.innerbin
outbin = snakemake.input.outerbin
nighres.laminar.volumetric_layering(inbin.replace('.nii.gz','p2l-surf.nii.gz'), outbin.replace('.nii.gz','p2l-surf.nii.gz'),, save_data=True, output_dir=path)
print(f'levelset to ', file=logfile, flush=True)

## copy to output (volumetric_layering produces multiple files so we need to specify the right one)
copyfile(path + 'innersurf_bin_p2l-surf_layering-depth.nii.gz', snakemake.output.coords)
