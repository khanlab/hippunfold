"""
QC of surface topology and mesh quality.

The hippocampal (and dentate) surface should be a compact, orientable
2-manifold with boundary, homeomorphic to a closed disk: one connected
component, genus g = 0 and b = 1 boundary loop, so the Euler characteristic

    chi = V - E + F = 2c - 2g - b = 1   (c = number of components)

Handles (g > 0), e.g. from grey matter bridges between adjacent folds, are
flagged, as are extra boundary loops (holes, b > 1), which are treated as free
boundaries by the Laplace-Beltrami solve. Mesh quality (self-intersections,
degenerate and sliver triangles) is checked on the midthickness, inner and
outer surfaces, and the unfolded (AP, PD) mapping is checked for flipped
triangles and areal distortion.
"""

import html
from collections import Counter

import numpy as np
from lib.qc import ERROR, WARNING, QCResult
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

QC_NAME = "surface"

# default thresholds, can be overridden by the rule params. Defect counts
# default to 0, so that every defect is reported (as a warning).
DEFAULT_THRESHOLDS = {
    "max_genus": 0,
    "max_holes": 0,  # boundary loops in addition to the outer one
    "max_nonmanifold": 0,  # non-manifold edges + vertices
    "max_misoriented_faces": 0,
    "max_degenerate_faces": 0,  # zero-area or duplicate, per surface
    "max_self_intersecting_faces": 0,  # per surface
    "sliver_angle": 5.0,  # degrees, min triangle angle below this is a sliver
    "max_sliver_faces": 0,  # per surface
    "max_boundary_chords": 0,  # interior edges joining two boundary vertices
    "max_flipped_faces": 0,  # in the unfolded mapping
    "max_collapsed_faces": 0,  # zero area in the unfolded mapping
    "distortion_log2": 2.0,  # |log2 areal ratio| above this is distorted (4x)
    "max_distorted_fraction": 0.05,
}

# cap on the number of element indices stored per defect, and on defects listed
MAX_STORED = 5000
MAX_DEFECTS = 500


def _cap(idx):
    return np.asarray(idx)[:MAX_STORED]


# --- topology ---


def unique_edges(faces):
    """
    Unique undirected edges, and for each face edge (F*3, in the order
    (v0,v1), (v1,v2), (v2,v0)) the index of its unique edge.
    """
    directed = np.stack([faces, np.roll(faces, -1, axis=1)], axis=2).reshape(-1, 2)
    edges, inverse = np.unique(np.sort(directed, axis=1), axis=0, return_inverse=True)
    return edges, inverse.ravel(), directed


def _components(n, i, j):
    graph = sparse.coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    return connected_components(graph, directed=False)


def boundary_loops(edges, edge_counts):
    """
    Boundary loops, as arrays of vertices (unordered), from the connected
    components of the boundary edges (edges used by one face).
    """
    bedges = edges[edge_counts == 1]
    if len(bedges) == 0:
        return []
    verts, local = np.unique(bedges, return_inverse=True)
    local = local.reshape(-1, 2)
    _, comp = _components(len(verts), local[:, 0], local[:, 1])
    return [verts[comp == k] for k in range(comp.max() + 1)]


def orientation(faces, edges, inverse, directed, edge_counts):
    """
    Whether the mesh is orientable, and the faces whose winding disagrees with
    the majority of their connected component (if orientable).

    Faces are linked across manifold edges. Each face has two copies, one per
    winding: across a consistently wound edge (opposite directions in the two
    faces) like copies are linked, otherwise opposite copies. The mesh is
    non-orientable if the two copies of a face are connected.
    """
    n_faces = len(faces)
    face_of = np.repeat(np.arange(n_faces), 3)
    manifold = edge_counts[inverse] == 2
    order = np.argsort(inverse[manifold], kind="stable")
    pairs = np.flatnonzero(manifold)[order].reshape(-1, 2)  # face-edge pairs
    f1, f2 = face_of[pairs[:, 0]], face_of[pairs[:, 1]]
    same_dir = np.all(directed[pairs[:, 0]] == directed[pairs[:, 1]], axis=1)
    # copy 0 of face f is node f, copy 1 is node f + n_faces
    i = np.r_[f1, f1 + n_faces]
    j = np.r_[f2 + n_faces * same_dir, f2 + n_faces * (1 - same_dir)]
    _, comp = _components(2 * n_faces, i, j)
    orientable = not np.any(comp[:n_faces] == comp[n_faces:])

    misoriented = np.array([], dtype=int)
    if orientable:
        # sign of each face relative to its component's reference winding
        sign = comp[:n_faces] < comp[n_faces:]
        _, face_comp = _components(n_faces, f1, f2)
        for k in np.unique(face_comp):
            members = np.flatnonzero(face_comp == k)
            s = sign[members]
            minority = s if s.sum() * 2 < len(s) else ~s
            misoriented = np.r_[misoriented, members[minority]]
    return orientable, int(np.sum(same_dir)), misoriented


def corner_groups(faces, inverse, edge_counts):
    """
    Group label of each face corner (corner c = 3*f + k is vertex faces[f, k]).
    Corners of a vertex are linked across manifold edges containing that
    vertex, so each group is one fan of faces around the vertex.
    """
    n_faces = len(faces)
    corners = faces.ravel()  # corner c = 3*f + k is vertex faces[f, k]
    # face edge e = 3*f + k joins corners 3*f + k and 3*f + (k+1)%3
    fe = np.arange(3 * n_faces)
    fe_c1 = fe
    fe_c2 = 3 * (fe // 3) + (fe % 3 + 1) % 3
    manifold = edge_counts[inverse] == 2
    order = np.argsort(inverse[manifold], kind="stable")
    pairs = np.flatnonzero(manifold)[order].reshape(-1, 2)
    a, b = pairs[:, 0], pairs[:, 1]
    # match corners on the same vertex between the two faces of each edge
    i, j = [], []
    for ca in (fe_c1[a], fe_c2[a]):
        for cb in (fe_c1[b], fe_c2[b]):
            same = corners[ca] == corners[cb]
            i.append(ca[same])
            j.append(cb[same])
    _, comp = _components(3 * n_faces, np.concatenate(i), np.concatenate(j))
    return comp


def nonmanifold_vertices(faces, inverse, edge_counts):
    """
    Vertices whose incident faces do not form a single fan (e.g. two sheets
    touching at a vertex), i.e. with more than one corner group.
    """
    comp = corner_groups(faces, inverse, edge_counts)
    groups = np.unique(np.c_[faces.ravel(), comp], axis=0)
    verts, n_groups = np.unique(groups[:, 0], return_counts=True)
    return verts[n_groups > 1]


def compute_topology(points, faces):
    """Topological invariants and defects of a triangle mesh."""
    faces = np.asarray(faces, dtype=np.int64)
    edges, inverse, directed = unique_edges(faces)
    edge_counts = np.bincount(inverse, minlength=len(edges))
    used = np.unique(faces)
    V, E, F = len(used), len(edges), len(faces)
    chi = V - E + F

    n_comp, _ = _components(len(points), edges[:, 0], edges[:, 1])
    n_comp -= len(points) - V  # unreferenced vertices are their own components

    loops = boundary_loops(edges, edge_counts)
    # largest loop (by perimeter) is the outer boundary, the rest are holes
    perimeters = []
    bedges = edges[edge_counts == 1]
    for loop in loops:
        in_loop = np.isin(bedges[:, 0], loop)
        seg = points[bedges[in_loop, 0]] - points[bedges[in_loop, 1]]
        perimeters.append(float(np.linalg.norm(seg, axis=1).sum()))
    order = np.argsort(perimeters)[::-1]
    loops = [loops[k] for k in order]
    perimeters = [perimeters[k] for k in order]

    nm_edges = edges[edge_counts > 2]
    nm_verts = nonmanifold_vertices(faces, inverse, edge_counts)
    orientable, n_inconsistent, misoriented = orientation(
        faces, edges, inverse, directed, edge_counts
    )
    b = len(loops)
    genus = (2 * n_comp - b - chi) / 2
    return {
        "V": V,
        "E": E,
        "F": F,
        "chi": chi,
        "components": n_comp,
        "boundary_loops": b,
        "genus": genus,
        "loops": loops,
        "loop_perimeters": perimeters,
        "nonmanifold_edges": nm_edges,
        "nonmanifold_vertices": nm_verts,
        "unreferenced_vertices": len(points) - V,
        "orientable": orientable,
        "inconsistent_edges": n_inconsistent,
        "misoriented_faces": misoriented,
    }


# --- triangle quality ---


def face_geometry(points, faces):
    """Face areas, minimum angles (degrees) and edge lengths (F x 3)."""
    tri = points[faces]
    e = np.stack(
        [tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 1], tri[:, 0] - tri[:, 2]], axis=1
    )
    lengths = np.linalg.norm(e, axis=2)
    area = 0.5 * np.linalg.norm(np.cross(e[:, 0], -e[:, 2]), axis=1)
    angles = []
    for k in range(3):
        u, v = -e[:, k - 1], e[:, k]  # edges leaving corner k
        cos = np.einsum("ij,ij->i", u, v) / (lengths[:, k - 1] * lengths[:, k] + 1e-20)
        angles.append(np.degrees(np.arccos(np.clip(cos, -1, 1))))
    return area, np.min(angles, axis=0), lengths


def max_angles(points, faces):
    """Largest angle (degrees) of each face."""
    tri = points[faces]
    angles = []
    for k in range(3):
        u = tri[:, (k + 1) % 3] - tri[:, k]
        v = tri[:, (k + 2) % 3] - tri[:, k]
        cos = np.einsum("ij,ij->i", u, v) / (
            np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) + 1e-20
        )
        angles.append(np.degrees(np.arccos(np.clip(cos, -1, 1))))
    return np.max(angles, axis=0)


def boundary_chords(faces, return_faces=False):
    """
    Interior edges (as vertex pairs) joining two boundary vertices, and
    optionally the faces they belong to.
    """
    faces = np.asarray(faces, dtype=np.int64)
    edges, inverse, _ = unique_edges(faces)
    counts = np.bincount(inverse, minlength=len(edges))
    boundary = np.zeros(faces.max() + 1, dtype=bool)
    boundary[edges[counts == 1].ravel()] = True
    chord = (counts == 2) & boundary[edges[:, 0]] & boundary[edges[:, 1]]
    if not return_faces:
        return edges[chord]
    return edges[chord], np.unique(np.flatnonzero(chord[inverse]) // 3)


def duplicate_faces(faces):
    _, counts = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    return int(np.sum(counts - 1))


# --- self-intersections ---


def _segments_hit_triangles(p0, p1, a, b, c, eps=1e-9):
    """Vectorised segment (p0->p1) vs triangle (a,b,c) intersection, strict."""
    d = p1 - p0
    e1, e2 = b - a, c - a
    h = np.cross(d, e2)
    det = np.einsum("ij,ij->i", e1, h)
    ok = np.abs(det) > eps
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = p0 - a
    u = inv * np.einsum("ij,ij->i", s, h)
    q = np.cross(s, e1)
    v = inv * np.einsum("ij,ij->i", d, q)
    t = inv * np.einsum("ij,ij->i", e2, q)
    return ok & (u > eps) & (v > eps) & (u + v < 1 - eps) & (t > eps) & (t < 1 - eps)


def self_intersecting_faces(points, faces, chunk=200000):
    """
    Faces that intersect another face they share no vertex with.
    Candidate pairs have overlapping bounding spheres (centroid, max distance
    to a corner); each pair is tested edge-vs-triangle in both directions.
    """
    points = np.asarray(points, dtype=float)
    tri = points[faces]
    centroid = tri.mean(axis=1)
    radius = np.linalg.norm(tri - centroid[:, None], axis=2).max(axis=1)
    tree = cKDTree(centroid)
    pairs = tree.query_pairs(2 * radius.max(), output_type="ndarray")
    if len(pairs) == 0:
        return np.array([], dtype=int)
    i, j = pairs[:, 0], pairs[:, 1]
    keep = np.linalg.norm(centroid[i] - centroid[j], axis=1) < radius[i] + radius[j]
    i, j = i[keep], j[keep]
    shared = (faces[i][:, :, None] == faces[j][:, None, :]).any(axis=(1, 2))
    i, j = i[~shared], j[~shared]

    hit = np.zeros(len(i), dtype=bool)
    for start in range(0, len(i), chunk):
        ii, jj = i[start : start + chunk], j[start : start + chunk]
        h = np.zeros(len(ii), dtype=bool)
        for x, y in ((ii, jj), (jj, ii)):
            for k in range(3):
                h |= _segments_hit_triangles(
                    tri[x, k], tri[x, (k + 1) % 3], tri[y, 0], tri[y, 1], tri[y, 2]
                )
        hit[start : start + chunk] = h
    return np.unique(np.r_[i[hit], j[hit]])


# --- unfolded mapping ---


def unfold_distortion(points, faces, coords_ap, coords_pd):
    """
    Flipped faces and areal distortion of the unfolded (AP, PD) mapping.

    Returns signed unfolded area, the collapsed faces (zero unfolded area,
    expected where all three vertices lie in a terminal band with a constant
    coordinate), the flipped faces (winding opposite to the majority) and log2
    of the ratio of each face's share of the unfolded area to its share of the
    native area (0 = no distortion, 1 = twice as large, -1 = half as large;
    NaN for collapsed faces).
    """
    uv = np.c_[coords_pd, coords_ap].astype(float)
    t = uv[faces]
    signed = 0.5 * (
        (t[:, 1, 0] - t[:, 0, 0]) * (t[:, 2, 1] - t[:, 0, 1])
        - (t[:, 2, 0] - t[:, 0, 0]) * (t[:, 1, 1] - t[:, 0, 1])
    )
    majority = np.sign(np.sum(np.sign(signed))) or 1.0
    tiny = np.abs(signed) <= 1e-12 * np.abs(signed).max()
    collapsed = np.flatnonzero(tiny)
    flipped = np.flatnonzero((np.sign(signed) != majority) & ~tiny)

    area3d, _, _ = face_geometry(np.asarray(points, dtype=float), faces)
    share2d = np.abs(signed) / np.abs(signed).sum()
    share3d = area3d / area3d.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        log_ratio = np.log2(share2d / share3d)
    log_ratio[tiny | ~np.isfinite(log_ratio)] = np.nan
    return signed, collapsed, flipped, log_ratio


# --- defects ---

# defect types: (category, description). Categories set the colour in figures.
DEFECT_TYPES = {
    "hole": ("topology", "Boundary loop in addition to the outer boundary"),
    "nonmanifold": ("topology", "Non-manifold edges or vertices"),
    "misoriented": ("topology", "Faces wound opposite to the rest of the surface"),
    "self_intersection": ("intersection", "Faces crossing other faces"),
    "degenerate": ("quality", "Zero-area or duplicate faces"),
    "sliver": ("quality", "Triangles with a very small angle"),
    "boundary_chord": (
        "quality",
        (
            "Interior edges joining two boundary vertices (faces with all "
            "vertices on one terminal arc collapse in the unfolded space)"
        ),
    ),
    "unfold_flipped": ("unfold", "Faces folded over in the unfolded mapping"),
    "unfold_collapsed": ("unfold", "Faces with zero area in the unfolded mapping"),
}
DEFECT_CATEGORIES = {"topology": 1, "intersection": 3, "quality": 2, "unfold": 0}


def face_clusters(faces, face_ids, n_vertices):
    """Group faces into clusters connected through shared vertices."""
    face_ids = np.unique(np.asarray(face_ids, dtype=np.int64))
    if len(face_ids) == 0:
        return []
    rows = np.repeat(np.arange(len(face_ids)), 3)
    incidence = sparse.csr_matrix(
        (np.ones(len(rows)), (rows, faces[face_ids].ravel())),
        shape=(len(face_ids), n_vertices),
    )
    n, comp = connected_components(incidence @ incidence.T, directed=False)
    clusters = [face_ids[comp == k] for k in range(n)]
    return sorted(clusters, key=len, reverse=True)


def _round_or_none(value, ndigits=3):
    return None if np.isnan(value) else round(float(value), ndigits)


def _vertex_faces(faces, vertices):
    return np.flatnonzero(np.isin(faces, vertices).any(axis=1))


class DefectList:
    """Defects found on a surface, each located on the midthickness surface."""

    def __init__(self, points, faces, coords_ap, coords_pd, boundary_vertices):
        """boundary_vertices: vertices of the outer boundary loop"""
        self.points = points
        self.faces = faces
        self.ap = np.asarray(coords_ap, dtype=float)
        self.pd = np.asarray(coords_pd, dtype=float)
        self.boundary = np.zeros(len(points), dtype=bool)
        self.boundary[boundary_vertices] = True
        self.defects = []

    def add(self, kind, surface, face_ids, detail="", vertices=None):
        """Add one defect for each connected cluster of the faces."""
        for cluster in face_clusters(self.faces, face_ids, len(self.points)):
            self.add_one(kind, surface, cluster, detail, vertices)

    def add_one(self, kind, surface, cluster, detail="", vertices=None):
        cluster = np.asarray(cluster, dtype=np.int64)
        verts = np.unique(self.faces[cluster]) if vertices is None else vertices
        centroids = self.points[self.faces[cluster]].mean(axis=1)
        # locate at the face nearest the cluster centroid, so it is on the surface
        at = cluster[np.argmin(np.linalg.norm(centroids - centroids.mean(0), axis=1))]
        fv = self.faces[at]
        d = detail(cluster) if callable(detail) else detail
        self.defects.append(
            {
                "type": kind,
                "category": DEFECT_TYPES[kind][0],
                "surface": surface,
                "n_faces": len(cluster),
                "n_vertices": len(verts),
                "location": np.round(self.points[fv].mean(axis=0), 2),
                "location_face": int(at),
                "ap": _round_or_none(self.ap[fv].mean()),
                "pd": _round_or_none(self.pd[fv].mean()),
                "at_boundary": bool(self.boundary[verts].any()),
                "detail": d,
                "faces": _cap(cluster),
                "vertices": _cap(verts),
            }
        )

    def to_list(self):
        order = list(DEFECT_TYPES)
        defects = sorted(
            self.defects, key=lambda d: (order.index(d["type"]), -d["n_faces"])
        )
        for i, d in enumerate(defects):
            d["id"] = i + 1
        return defects


# --- QC ---


def _fmt_list(values, n=5):
    values = list(values)
    s = ", ".join(f"{v:.1f}" for v in values[:n])
    return s + (", …" if len(values) > n else "")


def _per_surface(counts):
    return ", ".join(f"{name} {n}" for name, n in counts.items())


def _band(ap, pd):
    """Name the terminal band(s) where a coordinate is constant at 0 or 1."""
    names = []
    for coord, src, sink in ((pd, "PD src", "PD sink"), (ap, "AP src", "AP sink")):
        if np.ptp(coord) < 1e-9 and coord[0] in (0.0, 1.0):
            names.append(src if coord[0] == 0.0 else sink)
    return ", ".join(names)


def compute_qc(surfaces, coords_ap, coords_pd, entities, thresholds):
    """
    Metrics, checks and a list of defects for surface topology (midthickness),
    mesh quality and self-intersections (each surface) and the unfolded
    mapping.

    surfaces maps a surface name (midthickness first, then e.g. inner, outer,
    with the same vertices and faces) to (points, faces). Without coords (None),
    e.g. for a surface before the Laplace-Beltrami solve, the unfolded mapping
    is not checked. Returns a QCResult.
    """
    th = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    result = QCResult(QC_NAME, entities)
    names = list(surfaces)
    points, faces = surfaces[names[0]]
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    has_coords = coords_ap is not None and coords_pd is not None
    if has_coords:
        coords_ap = np.asarray(coords_ap, dtype=float)
        coords_pd = np.asarray(coords_pd, dtype=float)
    else:
        coords_ap = coords_pd = np.full(len(points), np.nan)

    # --- topology ---
    topo = compute_topology(points, faces)
    for key, desc in [
        ("V", "Vertices (referenced by a face)"),
        ("E", "Edges"),
        ("F", "Faces"),
        ("chi", "Euler characteristic V - E + F (1 for a disk)"),
        ("components", "Connected components c (1 for a disk)"),
        ("boundary_loops", "Boundary loops b (1 for a disk)"),
        ("genus", "Genus g = (2c - b - chi) / 2 (0 for a disk)"),
    ]:
        result.add_metric(key, topo[key], desc)
    result.add_metric(
        "boundary_perimeters",
        [round(p, 2) for p in topo["loop_perimeters"]],
        "Perimeter (mm) of each boundary loop, outer first",
    )
    n_nm_e = len(topo["nonmanifold_edges"])
    n_nm_v = len(topo["nonmanifold_vertices"])
    n_mis = len(topo["misoriented_faces"])
    result.add_metric("nonmanifold_edges", n_nm_e, "Edges with > 2 faces")
    result.add_metric(
        "nonmanifold_vertices", n_nm_v, "Vertices whose faces are not a single fan"
    )
    result.add_metric(
        "unreferenced_vertices",
        topo["unreferenced_vertices"],
        "Vertices not used by any face",
    )
    result.add_metric(
        "misoriented_faces", n_mis, "Faces wound opposite to the rest of the surface"
    )

    c, b, g, chi = (
        topo["components"],
        topo["boundary_loops"],
        topo["genus"],
        topo["chi"],
    )
    result.add_check(
        "single_component",
        c == 1,
        f"{c} connected components",
        severity=ERROR,
        description="The surface is one connected piece",
    )
    result.add_check(
        "has_boundary",
        b >= 1,
        f"{b} boundary loops" + (" (closed surface, no terminals)" if b == 0 else ""),
        severity=ERROR,
        description="The surface has a boundary, for the src/sink terminals",
    )
    result.add_check(
        "orientable",
        topo["orientable"] and n_mis <= th["max_misoriented_faces"],
        (
            "Non-orientable surface"
            if not topo["orientable"]
            else f"{n_mis} faces wound opposite to the rest"
        ),
        severity=ERROR if not topo["orientable"] else WARNING,
        description="The surface is orientable and its faces are consistently wound",
    )
    result.add_check(
        "manifold",
        n_nm_e + n_nm_v <= th["max_nonmanifold"],
        f"{n_nm_e} non-manifold edges, {n_nm_v} non-manifold vertices",
        severity=WARNING,
        description="Every edge has 1 or 2 faces and every vertex a single fan of "
        "faces. Otherwise genus and boundary loops are not well defined.",
    )
    result.add_check(
        "genus",
        g <= th["max_genus"],
        f"χ = {topo['V']} − {topo['E']} + {topo['F']} = {chi}; c = {c}, b = {b}, "
        f"g = {g:g}"
        + ("" if g <= th["max_genus"] else " — handle(s), e.g. a bridge between folds"),
        severity=WARNING,
        description="No handles (genus 0), as expected for a disk",
    )
    n_holes = max(b - 1, 0)
    result.add_check(
        "holes",
        n_holes <= th["max_holes"],
        f"{n_holes} holes in addition to the outer boundary"
        + (
            f", perimeters (mm): {_fmt_list(topo['loop_perimeters'][1:])}"
            if n_holes
            else ""
        ),
        severity=WARNING,
        description="One boundary loop. Holes are treated as free boundaries in "
        "the Laplace-Beltrami solve.",
    )

    loops = topo["loops"]
    defects = DefectList(
        points,
        faces,
        coords_ap,
        coords_pd,
        loops[0] if loops else np.array([], dtype=int),  # outer boundary
    )
    for loop, perimeter in zip(loops[1:], topo["loop_perimeters"][1:]):
        defects.add_one(
            "hole",
            names[0],
            _vertex_faces(faces, loop),
            f"perimeter {perimeter:.1f} mm",
            vertices=loop,
        )
    nm_verts = np.unique(
        np.r_[topo["nonmanifold_vertices"], topo["nonmanifold_edges"].ravel()]
    ).astype(np.int64)
    defects.add("nonmanifold", names[0], _vertex_faces(faces, nm_verts))
    defects.add(
        "misoriented",
        names[0],
        topo["misoriented_faces"],
        lambda cl: (
            "faces "
            + ", ".join(str(f) for f in cl[:6])
            + (", …" if len(cl) > 6 else "")
        ),
    )

    chords, chord_faces = boundary_chords(faces, return_faces=True)
    result.add_metric(
        "boundary_chords",
        len(chords),
        "Interior edges joining two boundary vertices",
    )
    result.add_check(
        "boundary_chords",
        len(chords) <= th["max_boundary_chords"],
        f"{len(chords)} interior edges join two boundary vertices",
        severity=WARNING,
        description="Every face has a vertex off the boundary, so none collapses "
        "when the terminal arcs of the boundary are fixed",
    )
    defects.add("boundary_chord", names[0], chord_faces)

    # --- mesh quality and self-intersections, per surface ---
    counts = {"degenerate": {}, "sliver": {}, "self_intersection": {}}
    min_angles = {}
    for name in names:
        p = np.asarray(surfaces[name][0], dtype=float)
        f = np.asarray(surfaces[name][1], dtype=np.int64)
        same_mesh = p.shape == points.shape and np.array_equal(f, faces)
        area, min_angle, lengths = face_geometry(p, f)
        min_angles[name] = min_angle
        suffix = "" if name == names[0] else f"_{name}"
        label = f"{name} surface"

        result.add_metric(
            f"surface_area{suffix}", float(area.sum()), f"Area (mm²) of the {label}"
        )
        result.add_metric(
            f"edge_length_median{suffix}",
            float(np.median(lengths)),
            f"Median edge length (mm) of the {label}",
        )
        result.add_metric(
            f"min_angle_p1{suffix}",
            float(np.percentile(min_angle, 1)),
            f"1st percentile of triangle min angle (°), {label}",
        )
        max_angle = max_angles(p, f)
        result.add_metric(
            f"obtuse_fraction{suffix}",
            float(np.mean(max_angle > 90)),
            f"Fraction of faces with an angle > 90° (negative cotangent weight), "
            f"{label}",
        )
        result.add_metric(
            f"max_angle_p99{suffix}",
            float(np.percentile(max_angle, 99)),
            f"99th percentile of triangle max angle (°), {label}",
        )

        repeated = (f[:, 0] == f[:, 1]) | (f[:, 1] == f[:, 2]) | (f[:, 0] == f[:, 2])
        degenerate = np.flatnonzero((area <= 1e-6 * np.median(area)) | repeated)
        n_dup = duplicate_faces(f)
        slivers = np.flatnonzero(min_angle < th["sliver_angle"])
        hits = self_intersecting_faces(p, f)
        counts["degenerate"][name] = len(degenerate) + n_dup
        counts["sliver"][name] = len(slivers)
        counts["self_intersection"][name] = len(hits)
        result.add_metric(
            f"degenerate_faces{suffix}", len(degenerate), f"Zero-area faces, {label}"
        )
        result.add_metric(
            f"duplicate_faces{suffix}", n_dup, f"Duplicate faces, {label}"
        )
        result.add_metric(
            f"sliver_faces{suffix}",
            len(slivers),
            f"Faces with min angle < {th['sliver_angle']:g}°, {label}",
        )
        result.add_metric(
            f"self_intersecting_faces{suffix}",
            len(hits),
            f"Faces crossing another face, {label}",
        )

        if same_mesh:
            defects.add("degenerate", name, degenerate)
            defects.add(
                "sliver",
                name,
                slivers,
                lambda cl, a=min_angle: f"min angle {a[cl].min():.1f}°",
            )
            defects.add("self_intersection", name, hits)

    result.add_check(
        "degenerate_faces",
        max(counts["degenerate"].values()) <= th["max_degenerate_faces"],
        f"Zero-area or duplicate faces: {_per_surface(counts['degenerate'])}",
        severity=WARNING,
        description="No zero-area or duplicated triangles",
    )
    result.add_check(
        "sliver_faces",
        max(counts["sliver"].values()) <= th["max_sliver_faces"],
        f"Faces with min angle < {th['sliver_angle']:g}°: "
        f"{_per_surface(counts['sliver'])}",
        severity=WARNING,
        description="No sliver triangles, which degrade the Laplace-Beltrami solve "
        "and curvature",
    )
    result.add_check(
        "self_intersections",
        max(counts["self_intersection"].values()) <= th["max_self_intersecting_faces"],
        f"Faces crossing another face: {_per_surface(counts['self_intersection'])}",
        severity=WARNING,
        description="No surface intersects itself, e.g. where adjacent folds touch",
    )

    # --- unfolded mapping ---
    if has_coords:
        _, collapsed, flipped, log_ratio = unfold_distortion(
            points, faces, coords_ap, coords_pd
        )
        finite = log_ratio[np.isfinite(log_ratio)]
        distorted = np.abs(finite) > th["distortion_log2"]
        result.add_metric(
            "unfold_collapsed_faces", len(collapsed), "Faces with zero unfolded area"
        )
        result.add_metric(
            "unfold_flipped_faces", len(flipped), "Flipped unfolded faces"
        )
        q25, q75 = np.percentile(finite, [25, 75])
        result.add_metric(
            "unfold_log2_distortion_iqr",
            float(q75 - q25),
            "IQR of log2 areal distortion, unfolded vs native",
        )
        result.add_metric(
            "unfold_distorted_fraction",
            float(distorted.mean()),
            f"Fraction of faces with |log2 areal distortion| > {th['distortion_log2']:g}",
        )
        in_band = sum(
            bool(_band(coords_ap[faces[i]], coords_pd[faces[i]])) for i in collapsed
        )
        result.add_check(
            "unfold_flipped",
            len(flipped) <= th["max_flipped_faces"],
            f"{len(flipped)} faces ({len(flipped) / len(faces):.2%}) folded over in the "
            "unfolded (AP, PD) mapping",
            severity=WARNING,
            description="The unfolded mapping does not fold over itself",
        )
        result.add_check(
            "unfold_collapsed",
            len(collapsed) <= th["max_collapsed_faces"],
            f"{len(collapsed)} faces ({len(collapsed) / len(faces):.2%}) have zero "
            f"unfolded area, {in_band} of them with all vertices in one terminal band",
            severity=WARNING,
            description="No faces collapse to zero area in the unfolded mapping (they "
            "do when all three vertices are in a terminal band, where a coordinate is "
            "fixed)",
        )
        result.add_check(
            "unfold_distortion",
            distorted.mean() <= th["max_distorted_fraction"],
            f"{distorted.mean():.1%} of faces change area share by more than "
            f"{2 ** th['distortion_log2']:g}× (max {th['max_distorted_fraction']:.0%})",
            severity=WARNING,
            description="Areal distortion between native and unfolded space is moderate",
        )
        misoriented = set(topo["misoriented_faces"].tolist())
        defects.add(
            "unfold_flipped",
            "unfold",
            flipped,
            lambda cl: (
                f"{sum(f in misoriented for f in cl.tolist())} of these faces are "
                "misoriented (winding, not mapping)"
                if misoriented.intersection(cl.tolist())
                else ""
            ),
        )
        defects.add(
            "unfold_collapsed",
            "unfold",
            collapsed,
            lambda cl: (
                _band(coords_ap[faces[cl]].ravel(), coords_pd[faces[cl]].ravel())
                and f"in {_band(coords_ap[faces[cl]].ravel(), coords_pd[faces[cl]].ravel())} band"
            ),
        )

    defect_list = defects.to_list()
    result.add_metric("defects", len(defect_list), "Defects (connected clusters) found")
    # topological defects come first, so only e.g. slivers are left out
    result.add_data("defects", defect_list[:MAX_DEFECTS])
    result.add_data("loops", [np.asarray(loop) for loop in loops])
    result.add_data("surfaces", names)
    return result


# --- figures, used by the report ---


def _mark(ax, xy, color, size, ring=True, labels=None):
    xy = np.asarray(xy, dtype=float).reshape(-1, 2)
    if len(xy) == 0:
        return
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        s=size,
        facecolors="none" if ring else color,
        edgecolors=color,
        linewidths=1.2,
        zorder=3,
    )
    for (x, y), text in zip(xy, labels or []):
        ax.annotate(
            text,
            (x, y),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
            fontweight="bold",
            color=color,
            zorder=4,
        )


def defect_legend():
    """(label, series index) pairs for the defect categories."""
    return [("outer boundary", None)] + [
        (cat, i) for cat, i in DEFECT_CATEGORIES.items()
    ]


def plot_defects(result, points, faces):
    """
    Native surface viewed face-on, with the defect clusters coloured by
    category and numbered as in the defect table, and the boundary loops.
    """
    from lib.qc_report import (
        INK,
        SERIES,
        face_normals_z,
        hex_rgb,
        render_faces,
        view_frame,
    )

    p = view_frame(np.asarray(points, dtype=float))
    faces = np.asarray(faces, dtype=np.int64)
    p2d = p[:, :2]
    face_rgb = np.tile([0.80, 0.80, 0.78], (len(faces), 1))
    series = hex_rgb(SERIES)
    defects = result.data.get("defects", [])
    for d in defects:
        face_rgb[np.asarray(d["faces"], dtype=int)] = series[
            DEFECT_CATEGORIES[d["category"]]
        ]
    fig = render_faces(
        p2d, p[:, 2], faces, face_rgb, face_normals_z(p, faces), figsize=(9, 6)
    )
    ax = fig.axes[0]

    loops = result.data.get("loops", [])
    if loops:
        _mark(ax, p2d[np.asarray(loops[0], dtype=int)], INK["text_secondary"], 1, False)
    centroid = p2d[faces].mean(axis=1)
    # with many defects (e.g. slivers of a raw isosurface), draw quality
    # defects as small dots without numbers, so the others stay readable
    crowded = len(defects) > 100
    for d in defects:
        color = SERIES[DEFECT_CATEGORIES[d["category"]]]
        at = centroid[[d["location_face"]]]
        if crowded and d["category"] == "quality":
            _mark(ax, at, color, 4, ring=False)
        else:
            _mark(ax, at, color, 45, labels=[str(d["id"])])
    return fig


def plot_distortion(points, faces, coords_ap, coords_pd, flipped, log_ratio, vmax=2.0):
    """
    Areal distortion (log2) on the native surface (face-on) and in unfolded
    (PD, AP) space, with flipped faces circled. Returns two figures.
    """
    from lib.qc_report import (
        DIVERGING_RAMP,
        INK,
        face_normals_z,
        hex_rgb,
        ramp_rgb,
        render_faces,
        view_frame,
    )

    faces = np.asarray(faces, dtype=np.int64)
    rgb = ramp_rgb(log_ratio, DIVERGING_RAMP, -vmax, vmax)
    rgb[~np.isfinite(log_ratio)] = hex_rgb([INK["muted"]])[0]  # collapsed
    flipped = np.asarray(flipped, dtype=int)

    p = view_frame(np.asarray(points, dtype=float))
    native = render_faces(p[:, :2], p[:, 2], faces, rgb, face_normals_z(p, faces))
    _mark(
        native.axes[0], p[:, :2][faces[flipped]].mean(axis=1), INK["text_primary"], 40
    )

    uv = np.c_[coords_pd, coords_ap].astype(float)
    unfolded = render_faces(uv, np.zeros(len(uv)), faces, rgb, np.ones(len(faces)))
    _mark(unfolded.axes[0], uv[faces[flipped]].mean(axis=1), INK["text_primary"], 40)
    return [native, unfolded]


def plot_histograms(min_angles, log_ratio, sliver_angle, distortion_log2):
    """
    Triangle min angle per surface (dict name -> angles), and log2 areal
    distortion of the unfolded mapping.
    """
    import matplotlib.pyplot as plt
    from lib.qc_report import INK, SERIES, style

    with style():
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.0))
        bins = np.linspace(0, 60, 61)
        for k, (name, angles) in enumerate(min_angles.items()):
            ax1.hist(
                angles,
                bins=bins,
                histtype="step",
                color=SERIES[k % len(SERIES)],
                lw=1.5,
                label=name,
                density=True,
            )
        ax1.axvline(sliver_angle, color=INK["muted"], lw=1, ls="--")
        ax1.set_xlabel("triangle min angle (°)")
        ax1.set_ylabel("density")
        ax1.legend(loc="upper left", fontsize=8)
        ax1.grid(axis="x", visible=False)

        finite = log_ratio[np.isfinite(log_ratio)]
        lim = max(3.0, distortion_log2 + 1)
        ax2.hist(
            np.clip(finite, -lim, lim),
            bins=np.linspace(-lim, lim, 61),
            color=SERIES[0],
            alpha=0.85,
        )
        for x in (-distortion_log2, distortion_log2):
            ax2.axvline(x, color=INK["muted"], lw=1, ls="--")
        ax2.set_xlabel("log2 areal distortion (unfolded / native)")
        ax2.set_ylabel("faces")
        ax2.grid(axis="x", visible=False)
        fig.subplots_adjust(wspace=0.25)
    return fig


# --- report helpers ---

# labels for the intermediate stages of processing, see surface_qc_stages
STAGE_LABELS = {"isosurf": "Isosurface", "repaired": "Repaired", "final": "Final"}
STAGE_DESCS = {
    "isosurf": "raw isosurface, before repair",
    "repaired": "after repair and decimation, before triangle-quality improvement",
}
# metrics compared between stages
COMPARE_METRICS = [
    "V",
    "F",
    "chi",
    "components",
    "boundary_loops",
    "genus",
    "nonmanifold_edges",
    "nonmanifold_vertices",
    "misoriented_faces",
    "unreferenced_vertices",
    "degenerate_faces",
    "duplicate_faces",
    "sliver_faces",
    "boundary_chords",
    "min_angle_p1",
    "obtuse_fraction",
    "max_angle_p99",
    "self_intersecting_faces",
    "surface_area",
    "edge_length_median",
    "defects",
]


def fmt_coords(d):
    if d["ap"] is None:
        return "–"
    return f"{d['ap']:.2f}, {d['pd']:.2f}"


def defect_table(defects, total=None):
    """HTML table of defects; total is the number found, if more are listed"""
    if not defects:
        return '<p class="muted">No defects found.</p>'
    rows = []
    for d in defects:
        series = DEFECT_CATEGORIES[d["category"]] + 1
        x, y, z = d["location"]
        rows.append(
            f'<tr><td class="num">{d["id"]}</td>'
            f'<td><span class="key"><span class="swatch" style="background:'
            f'var(--series-{series})"></span>{html.escape(d["type"])}</span></td>'
            f"<td>{html.escape(d['surface'])}</td>"
            f'<td class="num">{d["n_faces"]}</td>'
            f'<td class="num">{x:.1f}, {y:.1f}, {z:.1f}</td>'
            f'<td class="num">{fmt_coords(d)}</td>'
            f"<td>{'yes' if d['at_boundary'] else ''}</td>"
            f"<td>{html.escape(d['detail'] or '')}</td></tr>"
        )
    counts = Counter(d["type"] for d in defects)
    summary = ", ".join(f"{n} {t}" for t, n in counts.items())
    if total is not None and total > len(defects):
        summary += f" (first {len(defects)} of {total} listed)"
    return (
        f'<p class="muted" style="font-size:12px">{len(defects)} defects: {summary}.</p>'
        '<div class="table-wrap" style="max-height:440px;overflow:auto">'
        '<table><thead><tr><th class="num">#</th><th>Type</th><th>Surface</th>'
        '<th class="num">Faces</th><th class="num">Location (mm)</th>'
        '<th class="num">AP, PD</th><th>Outer boundary</th><th>Detail</th></tr></thead>'
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def defect_panel(title, result, points, faces):
    """Report panel with the defect map of a surface."""
    from lib.qc_report import fig_to_png, legend

    key = legend([(cat, i) for cat, i in DEFECT_CATEGORIES.items()])
    img = fig_to_png(plot_defects(result, points, faces))
    return (
        f'<div class="panel"><h4>{html.escape(title)}</h4>{key}'
        f'<img src="{img}" alt="Defect locations, {html.escape(title)}"></div>'
    )


DEFECT_MAP_CAPTION = (
    '<p class="muted" style="font-size:12px">Midthickness, face-on. Faces of each '
    "defect are coloured by category, and numbered as in the tables (defects on the "
    "inner or outer surface are shown at the same vertices of the midthickness). The "
    "outer boundary is outlined. Location is the centroid of the defect's central "
    "face, in the surface's space (corobl, mm); face and vertex ids of each defect "
    "are in the qcmetrics JSON.</p>"
)
