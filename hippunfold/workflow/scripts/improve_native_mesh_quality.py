import nibabel as nib
import numpy as np
from lib.qc_surface import compute_topology, face_geometry, self_intersecting_faces
from lib.surface import write_points_faces_to_gifti
from lib.surface_quality import improve_quality
from lib.utils import setup_logger

log_file = snakemake.log[0] if snakemake.log else None
logger = setup_logger(log_file)

TOPOLOGY_KEYS = ("V", "F", "chi", "components", "boundary_loops", "genus")


def summary(points, faces):
    t = compute_topology(points, faces)
    min_angle = face_geometry(points, faces)[1]
    return {
        **{k: t[k] for k in TOPOLOGY_KEYS},
        "misoriented": len(t["misoriented_faces"]),
        "nonmanifold": len(t["nonmanifold_edges"]) + len(t["nonmanifold_vertices"]),
        "self_intersecting": len(self_intersecting_faces(points, faces)),
        "slivers": int(np.sum(min_angle < 5)),
        "min_angle_p1": round(float(np.percentile(min_angle, 1)), 1),
    }


points, faces = nib.load(snakemake.input.surf_gii).agg_data()
points = np.asarray(points, dtype=float)
faces = np.asarray(faces, dtype=np.int64)
before = summary(points, faces)
logger.info(f"before: {before}")

new_points, new_faces, info = improve_quality(
    points,
    faces,
    rounds=snakemake.params.rounds,
    relax_iterations=snakemake.params.relax_iterations,
    step=snakemake.params.relax_step,
    flip_tolerance=snakemake.params.flip_tolerance,
)
after = summary(new_points, new_faces)
logger.info(f"after: {after}")
logger.info(
    f"{info['edge_flips']} edge flips; vertices moved (along the surface) "
    f"max {info['max_vertex_move']:.3f} mm, mean {info['mean_vertex_move']:.3f} mm; "
    f"distance between the surfaces max {info['max_surface_distance']:.3f} mm, "
    f"mean {info['mean_surface_distance']:.4f} mm"
)

# keep the input if anything got worse that this step must not change
problems = [
    k for k in (*TOPOLOGY_KEYS, "misoriented", "nonmanifold") if after[k] != before[k]
]
if after["self_intersecting"] > before["self_intersecting"]:
    problems.append("self_intersecting")
if problems:
    logger.warning(
        f"triangle-quality improvement changed {problems}; keeping the input surface"
    )
    new_points, new_faces = points, faces

write_points_faces_to_gifti(new_points, new_faces, snakemake.output.surf_gii)
