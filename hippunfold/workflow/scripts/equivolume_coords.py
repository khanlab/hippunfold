import os
import nighres
from shutil import copyfile
import sys

tmpdir = sys.argv[1]
innerbin = sys.argv[2]
outerbin = sys.argv[3]
output_coords = sys.argv[4]

print("start")

nighres_args = {"save_data": True, "output_dir": tmpdir, "overwrite": True}

## convert binarized edges to levelset surfaces
levelset_inner = nighres.surface.probability_to_levelset(innerbin, **nighres_args)
levelset_outer = nighres.surface.probability_to_levelset(outerbin, **nighres_args)

print("binarized files to levelsets complete")

## get isovolume solution
isovolume = nighres.laminar.volumetric_layering(
    levelset_inner["result"], levelset_outer["result"], **nighres_args
)

print("levelset to isovolume complete")

## copy to output (volumetric_layering produces multiple files so we need to specify the right one)
copyfile(isovolume["depth"], output_coords)
