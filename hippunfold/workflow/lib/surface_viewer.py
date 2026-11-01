"""
Self-contained interactive HTML viewer for hippunfold surfaces.

Shows each surface in anatomical (corobl) and unfolded space side by side,
coloured by any of its shape metrics or by subfield labels, with linked hover.
The page draws with WebGL directly (no external libraries or assets), so it
opens offline. With several subjects it adds a group overview: every subject's
unfolded map of the chosen metric on shared colour limits, for spotting
outliers. Typical use:

    viewer = SurfaceViewer("Surfaces", entities={"subject": "001"})
    viewer.add("001", "L", "hipp", native_gii, unfold_gii,
               metrics={"thickness": path, ...}, subfields=label_gii)
    viewer.save("surfaces.html")

Surfaces resampled to a template density share their triangles and unfolded
coordinates across subjects; these are stored once per page.
"""

import base64
import datetime
import hashlib
import html
import json
from pathlib import Path

import matplotlib
import nibabel as nib
import numpy as np
from lib.qc_report import _CSS

# display settings per metric; others use the defaults
METRIC_STYLE = {
    "gyrification": {
        "name": "gyrification",
        "cmap": "viridis",
        "desc": "Smoothing-based gyrification index: local area of the "
        "midthickness over that of a copy smoothed along the surface. "
        "1 is flat; folds narrower than the smoothing scale raise it.",
    },
    "legacygyrification": {
        "name": "legacy gyrification",
        "cmap": "viridis",
        "desc": "Ratio of native to unfolded surface area (the scaling factor "
        "when unfolding), as used by unfolded registration with older atlases.",
    },
    "thickness": {
        "name": "thickness",
        "cmap": "viridis",
        "unit": "mm",
        "desc": "Distance between the inner and outer surfaces (mm).",
    },
    "curvature": {
        "name": "curvature",
        "cmap": "RdBu_r",
        "symmetric": True,
        "desc": "Mean curvature of the smoothed midthickness.",
    },
    "myelin": {
        "name": "myelin",
        "cmap": "magma",
        "desc": "T1w/T2w ratio sampled between the inner and outer surfaces.",
    },
}
DEFAULT_STYLE = {"cmap": "viridis", "desc": ""}
COLORMAPS = ("viridis", "magma", "RdBu_r")


def _b64(array):
    return base64.b64encode(np.ascontiguousarray(array).tobytes()).decode()


def _read_surface(path):
    gii = nib.load(path)
    points = gii.agg_data("NIFTI_INTENT_POINTSET").astype(np.float64)
    faces = gii.agg_data("NIFTI_INTENT_TRIANGLE").astype(np.int64)
    return points, faces


def _quantize(values, lo, hi):
    """uint8 codes 0-254 over [lo, hi]; 255 marks missing values"""
    values = np.asarray(values, float)
    span = hi - lo if hi > lo else 1.0
    q = np.clip(np.round((values - lo) / span * 254), 0, 254)
    q[~np.isfinite(values)] = 255
    return q.astype(np.uint8)


def _limits(values, symmetric=False):
    """robust colour limits: 2nd-98th percentile (symmetric about 0 if asked)"""
    finite = np.asarray(values, float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return [0.0, 1.0]
    lo, hi = np.percentile(finite, [2, 98])
    if symmetric:
        a = max(abs(lo), abs(hi)) or 1.0
        lo, hi = -a, a
    if hi <= lo:
        hi = lo + 1.0
    return [float(lo), float(hi)]


def _colormap_luts(n=64):
    return {
        name: [
            matplotlib.colors.to_hex(c)[1:]
            for c in matplotlib.colormaps[name](np.linspace(0, 1, n))
        ]
        for name in COLORMAPS
    }


class SurfaceViewer:
    def __init__(self, title, entities=None, description=""):
        self.title = title
        self.entities = dict(entities or {})
        self.description = description
        self.surfaces = []

    def add(
        self, subject, hemi, label, native_surf, unfold_surf, metrics, subfields=None
    ):
        """
        Add one surface. native_surf and unfold_surf are GIFTI surfaces with the
        same vertices; metrics maps metric names to .shape.gii files; subfields
        is an optional .label.gii.
        """
        points, faces = _read_surface(native_surf)
        flat, flat_faces = _read_surface(unfold_surf)
        if len(flat) != len(points) or not np.array_equal(faces, flat_faces):
            raise ValueError(f"{native_surf} and {unfold_surf} differ in topology")
        flat = flat[:, :2].copy()
        if hemi == "R":
            flat[:, 0] *= -1  # unfolded right hemispheres are mirrored in x
        values = {}
        for name, path in metrics.items():
            data = np.asarray(nib.load(path).agg_data(), float)
            if data.shape != (len(points),):
                raise ValueError(f"{path} does not match {native_surf}")
            values[name] = data
        labels = None
        if subfields is not None:
            gii = nib.load(subfields)
            keys = np.asarray(gii.agg_data(), int)
            table = {
                int(lab.key): [
                    lab.label,
                    matplotlib.colors.to_hex(lab.rgba[:3]),
                ]
                for lab in gii.labeltable.labels
                if lab.key in set(keys.tolist())
            }
            labels = (keys, table)
        self.surfaces.append(
            {
                "subject": str(subject),
                "hemi": hemi,
                "label": label,
                "points": points,
                "flat": flat,
                "faces": faces,
                "metrics": values,
                "labels": labels,
            }
        )

    def _data(self):
        blobs = {}

        def blob(array):
            raw = np.ascontiguousarray(array)
            key = hashlib.sha1(raw.tobytes() + str(raw.dtype).encode()).hexdigest()[:12]
            blobs.setdefault(key, _b64(raw))
            return key

        # shared colour limits per (label, metric), pooled over all surfaces
        shared = {}
        for s in self.surfaces:
            for name, v in s["metrics"].items():
                shared.setdefault((s["label"], name), []).append(v)
        shared = {
            f"{label}|{name}": _limits(
                np.concatenate(vs),
                METRIC_STYLE.get(name, DEFAULT_STYLE).get("symmetric", False),
            )
            for (label, name), vs in shared.items()
        }

        records = []
        for s in self.surfaces:
            points, faces = s["points"], s["faces"]
            lo, hi = points.min(0), points.max(0)
            span = np.where(hi > lo, hi - lo, 1.0)
            qpoints = np.round((points - lo) / span * 65535).astype(np.uint16)
            u32 = len(points) > 65535
            rec = {
                "subject": s["subject"],
                "hemi": s["hemi"],
                "label": s["label"],
                "nv": len(points),
                "nf": len(faces),
                "u32": bool(u32),
                "points": blob(qpoints),
                "origin": lo.tolist(),
                "span": span.tolist(),
                "flat": blob(s["flat"].astype(np.float32)),
                "faces": blob(faces.astype(np.uint32 if u32 else np.uint16)),
                "metrics": {},
            }
            for name, v in s["metrics"].items():
                style = METRIC_STYLE.get(name, DEFAULT_STYLE)
                finite = v[np.isfinite(v)]
                qlo, qhi = (
                    (float(finite.min()), float(finite.max()))
                    if finite.size
                    else (0.0, 1.0)
                )
                rec["metrics"][name] = {
                    "q": blob(_quantize(v, qlo, qhi)),
                    "qlo": qlo,
                    "qhi": qhi,
                    "own": _limits(v, style.get("symmetric", False)),
                    "median": float(np.median(finite)) if finite.size else None,
                }
            if s["labels"] is not None:
                keys, table = s["labels"]
                rec["subfields"] = {
                    "q": blob(np.clip(keys, 0, 255).astype(np.uint8)),
                    "table": table,
                }
            records.append(rec)

        metrics = []
        for s in self.surfaces:
            for name in s["metrics"]:
                if name not in metrics:
                    metrics.append(name)
        styles = {
            name: {**DEFAULT_STYLE, "name": name, **METRIC_STYLE.get(name, {})}
            for name in metrics
        }
        return {
            "records": records,
            "blobs": blobs,
            "shared": shared,
            "styles": styles,
            "metrics": metrics,
            "cmaps": _colormap_luts(),
        }

    def render(self):
        ents = " · ".join(
            f"{k}-{html.escape(str(v))}" for k, v in self.entities.items()
        )
        now = datetime.datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")
        n_subjects = len({s["subject"] for s in self.surfaces})
        desc = (
            f'<p class="lede">{html.escape(self.description)}</p>'
            if self.description
            else ""
        )
        title = self.title
        if self.entities:
            title += " — " + " ".join(f"{k}-{v}" for k, v in self.entities.items())
        page = _TEMPLATE.read_text()
        data = json.dumps(self._data(), separators=(",", ":")).replace("</", "<\\/")
        return (
            page.replace("{{css}}", _CSS)
            .replace("{{title}}", html.escape(title))
            .replace("{{heading}}", html.escape(self.title))
            .replace("{{entities}}", f"{ents} · " if ents else "")
            .replace("{{n_subjects}}", str(n_subjects))
            .replace("{{generated}}", now)
            .replace("{{description}}", desc)
            .replace("{{data}}", data)
        )

    def save(self, path):
        Path(path).write_text(self.render())


_TEMPLATE = Path(__file__).with_name("surface_viewer.html")
