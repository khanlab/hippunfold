import os
import nighres
from shutil import copyfile


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    print('start')

    nighres_args = {'save_data': True, 'output_dir': tmpdir, 'overwrite': True}

    tmpdir = snakemake.resources.tmpdir

    ## convert binarized edges to levelset surfaces
    levelset_inner = nighres.surface.probability_to_levelset(snakemake.input.innerbin, **nighres_args)
    levelset_outer = nighres.surface.probability_to_levelset(snakemake.input.outerbin, **nighres_args)

    print('binarized files to levelsets complete')


    ## get isovolume solution
    isovolume = nighres.laminar.volumetric_layering(levelset_inner['result'], levelset_outer['result'], **nighres_args)

    print('levelset to isovolume complete')


    ## copy to output (volumetric_layering produces multiple files so we need to specify the right one)
    copyfile(isovolume['depth'], snakemake.output.coords)

