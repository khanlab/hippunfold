import nibabel as nib
import numpy as np
from lib.gyrification import smoothing_gyrification
from lib.surface import read_surface_from_gifti
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
faces = surface.faces.reshape(-1, 4)[:, 1:]
scale = float(snakemake.params.scale)
sigma = float(snakemake.params.sigma)

logger.info(f"Smoothing-based gyrification: scale {scale} mm, pooling sigma {sigma} mm")
sgi = smoothing_gyrification(surface.points, faces, scale=scale, sigma=sigma)
logger.info(
    f"sGI median {np.median(sgi):.3f}, range {sgi.min():.3f} to {sgi.max():.3f}"
)

image = nib.GiftiImage()
image.meta["AnatomicalStructurePrimary"] = metadata["AnatomicalStructurePrimary"]
# record the scales, since they change the values
image.meta["GyrificationScale"] = str(scale)
image.meta["GyrificationPoolingSigma"] = str(sigma)
image.add_gifti_data_array(nib.gifti.GiftiDataArray(data=sgi.astype(np.float32)))
nib.save(image, snakemake.output.gii)
