"""QC of the AP/PD src/sink terminals assigned on the surface boundary."""

import numpy as np
from lib.qc import ERROR, WARNING, QCResult
from lib.surface import EXPECTED_TERMINAL_ORDER, TERMINAL_NAMES

QC_NAME = "srcsink"

# default thresholds, can be overridden by the rule params
DEFAULT_THRESHOLDS = {
    "max_relabelled_fraction": 0.1,
    "max_ap_fraction": 0.25,
    "max_terminal_distance": 3.0,
}


def _order_str(order):
    return " → ".join(TERMINAL_NAMES[k] for k in order)


def cyclic_runs(labels):
    """Number of runs of equal labels around a cyclic sequence."""
    labels = np.asarray(labels)
    n_changes = int(np.sum(labels != np.roll(labels, 1)))
    return max(n_changes, 1)


def compute_qc(distances, boundary_vertices, labels, info, hemi, entities, thresholds):
    """
    Metrics and checks for the terminals from assign_boundary_terminals.

    Returns a QCResult, with the loop data needed for plotting.
    """
    thresholds = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    result = QCResult(QC_NAME, entities)

    vertex_to_idx = {v: i for i, v in enumerate(boundary_vertices)}
    loop_idx = np.array([vertex_to_idx[v] for v in info["loop"]])
    loop_labels = np.asarray(info["loop_labels"])
    nearest = np.asarray(info["loop_nearest_labels"])
    loop_dist = np.asarray(distances)[loop_idx]
    m = len(loop_idx)

    result.add_metric("loop_vertices", m, "Vertices on the outer boundary loop")
    result.add_metric(
        "band_vertices", len(boundary_vertices), "Vertices in the boundary band"
    )
    result.add_metric(
        "terminal_order", _order_str(info["order"]), "Fitted order along the loop"
    )
    loop_counts = [int(np.sum(loop_labels == k)) for k in range(4)]
    for k, name in enumerate(TERMINAL_NAMES):
        result.add_metric(f"loop_vertices_{name}", loop_counts[k], f"{name} arc length")
    band_counts = [int(np.sum(np.asarray(labels) == k)) for k in range(4)]
    for k, name in enumerate(TERMINAL_NAMES):
        result.add_metric(
            f"band_vertices_{name}", band_counts[k], f"{name} band vertices"
        )

    # hemisphere-specific order
    expected = EXPECTED_TERMINAL_ORDER.get(hemi)
    if expected is None:
        result.add_check(
            "terminal_order_hemisphere",
            False,
            f"No expected terminal order for hemi-{hemi}, not checked",
            severity=WARNING,
            description="Terminal order along the loop matches the hemisphere",
        )
    else:
        result.add_check(
            "terminal_order_hemisphere",
            tuple(info["order"]) == tuple(expected),
            f"Order is {_order_str(info['order'])}; expected {_order_str(expected)} "
            f"for hemi-{hemi}",
            severity=ERROR,
            description="Terminal order along the loop (face winding) matches the "
            "hemisphere. A mismatch means swapped terminals or a flipped surface.",
        )

    # fragmentation of the nearest-terminal labels, before fitting
    n_runs = cyclic_runs(nearest)
    result.add_metric(
        "nearest_terminal_runs", n_runs, "Runs of nearest-terminal labels on the loop"
    )
    result.add_check(
        "nearest_terminal_fragmentation",
        n_runs == 4,
        f"{n_runs} runs of nearest-terminal labels along the loop (4 if clean"
        + (", fewer means a terminal is missing)" if n_runs < 4 else ")"),
        severity=WARNING,
        description="Labels from the volumetric terminals are already 4 clean arcs",
    )

    n_relabelled = int(np.sum(loop_labels != nearest))
    frac = n_relabelled / m
    result.add_metric(
        "relabelled_vertices", n_relabelled, "Loop vertices changed by the fit"
    )
    result.add_metric(
        "relabelled_fraction", frac, "Fraction of loop changed by the fit"
    )
    result.add_check(
        "relabelled_fraction",
        frac <= thresholds["max_relabelled_fraction"],
        f"{n_relabelled} of {m} loop vertices ({frac:.1%}) changed from the nearest "
        f"terminal (max {thresholds['max_relabelled_fraction']:.0%})",
        severity=WARNING,
        description="The fit only needed small changes to the nearest-terminal labels",
    )

    ap_fracs = [loop_counts[k] / m for k in (0, 1)]
    result.add_check(
        "ap_terminal_extent",
        max(ap_fracs) <= thresholds["max_ap_fraction"],
        f"AP src {ap_fracs[0]:.1%}, AP sink {ap_fracs[1]:.1%} of the loop "
        f"(max {thresholds['max_ap_fraction']:.0%})",
        severity=WARNING,
        description="AP terminals are short, point-like arcs",
    )

    medians = []
    for k in range(4):
        d = loop_dist[loop_labels == k, k]
        medians.append(float(np.median(d)))
        result.add_metric(
            f"median_distance_{TERMINAL_NAMES[k]}",
            medians[-1],
            f"Median distance (mm) of the {TERMINAL_NAMES[k]} arc to its terminal",
        )
    far = [
        f"{TERMINAL_NAMES[k]} {medians[k]:.2f} mm"
        for k in range(4)
        if medians[k] > thresholds["max_terminal_distance"]
    ]
    result.add_check(
        "terminal_distance",
        not far,
        (
            "Too far: " + ", ".join(far)
            if far
            else f"All arcs within {thresholds['max_terminal_distance']} mm "
            "(median) of their terminal"
        ),
        severity=WARNING,
        description="Each arc lies close to its volumetric terminal label",
    )

    result.add_check(
        "band_connected",
        info["n_unreached"] == 0,
        f"{info['n_unreached']} band vertices not connected to the loop",
        severity=WARNING,
        description="Every boundary band vertex is connected to the boundary loop",
    )

    result.add_data("loop", info["loop"])
    result.add_data("loop_labels", loop_labels)
    result.add_data("loop_nearest_labels", nearest)
    result.add_data("loop_distances", np.round(loop_dist, 3))
    result.add_data("boundary_vertices", boundary_vertices)
    result.add_data("labels", labels)
    result.add_data("order", info["order"])
    result.add_data("expected_order", expected)
    return result


# --- figures, used by the report ---


def plot_loop(result):
    """Distances along the loop, with fitted and nearest-terminal labels."""
    import matplotlib.pyplot as plt
    from lib.qc_report import SERIES, style

    dist = np.asarray(result.data["loop_distances"])
    fitted = np.asarray(result.data["loop_labels"])
    nearest = np.asarray(result.data["loop_nearest_labels"])
    m = len(fitted)
    x = np.arange(m)

    def runs(labels):
        starts = np.flatnonzero(np.r_[True, labels[1:] != labels[:-1]])
        ends = np.r_[starts[1:], m]
        return [(s, e - s, labels[s]) for s, e in zip(starts, ends)]

    with style():
        fig, (ax, ax_strip) = plt.subplots(
            2,
            1,
            figsize=(10, 3.8),
            sharex=True,
            gridspec_kw={"height_ratios": [4, 1], "hspace": 0.08},
        )
        # terminals are decided close to 0 mm, so don't let far values dominate
        ymax = min(10.0, max(1.0, float(np.percentile(dist, 90))))
        for s, w, k in runs(fitted):
            ax.axvspan(s - 0.5, s + w - 0.5, color=SERIES[k], alpha=0.10, lw=0)
        for k in range(4):
            ax.plot(x, dist[:, k], color=SERIES[k], lw=1.5, label=TERMINAL_NAMES[k])
        ax.set_ylim(min(0.0, float(dist.min())), ymax)
        ax.set_ylabel("distance to terminal (mm)")
        ax.grid(axis="x", visible=False)

        for row, labels in enumerate([fitted, nearest]):
            for s, w, k in runs(labels):
                ax_strip.broken_barh(
                    [(s - 0.5, w)], (1 - row - 0.4, 0.8), color=SERIES[k], lw=0
                )
        ax_strip.set_yticks([1, 0], ["fitted", "nearest"])
        ax_strip.set_ylim(-0.6, 1.6)
        ax_strip.grid(False)
        ax_strip.spines["left"].set_visible(False)
        ax_strip.tick_params(axis="y", length=0)
        ax_strip.set_xlim(-0.5, m - 0.5)
        ax_strip.set_xlabel(
            "position along boundary loop (vertices, face-winding order)"
        )
    return fig


def _view_frame(points):
    """Rotate points to their principal axes, viewing along the smallest."""
    centred = points - points.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    return centred @ vt.T


def _render_faces(points2d, depth, faces, face_rgb, normals_z):
    """Draw faces back to front with simple shading, no axes."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    shade = 0.55 + 0.45 * np.abs(normals_z)
    colors = np.clip(face_rgb * shade[:, None], 0, 1)
    order = np.argsort(depth[faces].mean(axis=1))
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.add_collection(
        PolyCollection(
            points2d[faces[order]],
            facecolors=colors[order],
            edgecolors="none",
            antialiased=False,
        )
    )
    ax.set_xlim(points2d[:, 0].min(), points2d[:, 0].max())
    ax.set_ylim(points2d[:, 1].min(), points2d[:, 1].max())
    ax.set_aspect("equal")
    ax.axis("off")
    return fig


def _hex_rgb(colors):
    return np.array([[int(c[i : i + 2], 16) / 255 for i in (1, 3, 5)] for c in colors])


def plot_surfaces(result, points, faces, coords_ap, coords_pd):
    """
    Figures of the surface (viewed face-on): terminal labels on the boundary,
    and the AP and PD Laplace coordinates. Returns a list of figures.
    """
    from lib.qc_report import BLUE_RAMP, ORANGE_RAMP, SERIES

    p = _view_frame(np.asarray(points, dtype=float))
    faces = np.asarray(faces)
    tri = p[faces]
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-12
    p2d, depth = p[:, :2], p[:, 2]

    # terminals: faces touching band vertices take the majority band label
    vertex_labels = np.full(len(p), -1)
    vertex_labels[np.asarray(result.data["boundary_vertices"])] = result.data["labels"]
    base = np.array([0.80, 0.80, 0.78])
    series = _hex_rgb(SERIES)
    face_labels = vertex_labels[faces]
    face_rgb = np.tile(base, (len(faces), 1))
    in_band = (face_labels >= 0).any(axis=1)
    majority = np.array(
        [np.bincount(fl[fl >= 0], minlength=4).argmax() for fl in face_labels[in_band]]
    )
    face_rgb[in_band] = series[majority]
    figs = [_render_faces(p2d, depth, faces, face_rgb, normals[:, 2])]

    for coords, ramp in ((coords_ap, BLUE_RAMP), (coords_pd, ORANGE_RAMP)):
        ramp_rgb = _hex_rgb(ramp)
        value = np.clip(np.asarray(coords)[faces].mean(axis=1), 0, 1)
        pos = value * (len(ramp) - 1)
        lo = np.floor(pos).astype(int).clip(0, len(ramp) - 2)
        t = (pos - lo)[:, None]
        rgb = ramp_rgb[lo] * (1 - t) + ramp_rgb[lo + 1] * t
        figs.append(_render_faces(p2d, depth, faces, rgb, normals[:, 2]))
    return figs
