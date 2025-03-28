from lib.surface import (
    read_surface_from_gifti,
    read_metric_from_gii,
    write_surface_to_gifti,
)


surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
ap = read_metric_from_gii(snakemake.input.coords_AP)
pd = read_metric_from_gii(snakemake.input.coords_AP)


surface.points[:, 0] = ap * -float(snakemake.params.vertspace["extent"][0]) - float(
    snakemake.params.vertspace["origin"][0]
)

surface.points[:, 1] = pd * -float(snakemake.params.vertspace["extent"][1]) - float(
    snakemake.params.vertspace["origin"][1]
)

surface.points[:, 2] = snakemake.params.z_level

# TODO: write metadata too

write_surface_to_gifti(surface, snakemake.output.surf_gii)
