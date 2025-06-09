from lib.surface import read_surface_from_gifti, write_surface_to_gifti

surface, metadata = read_surface_from_gifti(snakemake.input.surf_gii)
surface.points[:, 2] = snakemake.params.z_level
write_surface_to_gifti(surface, snakemake.output.surf_gii)
