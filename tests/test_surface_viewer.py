import base64
import json
import re
import shutil
import subprocess

import nibabel as nib
import numpy as np
import pytest
from lib.surface_viewer import SurfaceViewer
from test_gyrification import structured_sheet


def write_surface(path, points, faces):
    gii = nib.GiftiImage()
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            points.astype(np.float32), intent="NIFTI_INTENT_POINTSET"
        )
    )
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(faces.astype(np.int32), intent="NIFTI_INTENT_TRIANGLE")
    )
    nib.save(gii, path)


def write_metric(path, values):
    gii = nib.GiftiImage()
    gii.add_gifti_data_array(nib.gifti.GiftiDataArray(values.astype(np.float32)))
    nib.save(gii, path)


def write_labels(path, keys):
    gii = nib.GiftiImage()
    for key, name, rgba in [
        (0, "???", (1, 1, 1, 0)),
        (1, "CA1", (1, 0, 0, 1)),
        (2, "CA3", (0, 0, 1, 1)),
    ]:
        lab = nib.gifti.GiftiLabel(key, *rgba)
        lab.label = name
        gii.labeltable.labels.append(lab)
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(keys.astype(np.int32), intent="NIFTI_INTENT_LABEL")
    )
    nib.save(gii, path)


@pytest.fixture
def surfaces(tmp_path):
    """two subjects sharing the triangles and unfolded coordinates, as at den-8k"""
    xy, faces = structured_sheet(21)
    unfold = tmp_path / "unfold.surf.gii"
    write_surface(unfold, np.c_[xy, np.zeros(len(xy))], faces)
    out = []
    for k, sub in enumerate(["sub-01", "sub-02"]):
        native = tmp_path / f"{sub}_native.surf.gii"
        write_surface(native, np.c_[xy, np.sin(xy[:, 0] + k)], faces)
        thick = tmp_path / f"{sub}_thickness.shape.gii"
        write_metric(thick, 1.0 + k + 0.1 * xy[:, 1])
        labels = tmp_path / f"{sub}_subfields.label.gii"
        write_labels(labels, (xy[:, 0] > 6).astype(int) + 1)
        out.append((sub, native, unfold, {"thickness": thick}, labels))
    return out


def page_data(page):
    match = re.search(
        r'<script id="data" type="application/json">(.*?)</script>', page, re.DOTALL
    )
    return json.loads(match.group(1))


def decode(data, key, dtype):
    return np.frombuffer(base64.b64decode(data["blobs"][key]), dtype=dtype)


def test_page_is_self_contained(surfaces):
    viewer = SurfaceViewer("Surfaces", entities={"den": "8k"})
    for sub, native, unfold, metrics, labels in surfaces:
        viewer.add(sub, "L", "hipp", native, unfold, metrics, subfields=labels)
    page = viewer.render()
    assert "{{" not in page
    assert not re.search(
        r'<(script|link)[^>]+(src|href)="https?:', page
    )  # no external assets
    data = page_data(page)
    assert [r["subject"] for r in data["records"]] == ["sub-01", "sub-02"]


def test_shared_geometry_stored_once(surfaces):
    viewer = SurfaceViewer("Surfaces")
    for sub, native, unfold, metrics, labels in surfaces:
        viewer.add(sub, "L", "hipp", native, unfold, metrics, subfields=labels)
    a, b = page_data(viewer.render())["records"]
    assert a["faces"] == b["faces"] and a["flat"] == b["flat"]
    assert a["points"] != b["points"]


def test_values_round_trip(surfaces):
    sub, native, unfold, metrics, labels = surfaces[0]
    viewer = SurfaceViewer("Surfaces")
    viewer.add(sub, "R", "hipp", native, unfold, metrics, subfields=labels)
    data = page_data(viewer.render())
    rec = data["records"][0]

    points = nib.load(native).agg_data("NIFTI_INTENT_POINTSET")
    q = decode(data, rec["points"], np.uint16).reshape(-1, 3)
    decoded = np.array(rec["origin"]) + q / 65535 * np.array(rec["span"])
    np.testing.assert_allclose(decoded, points, atol=1e-3)

    flat = decode(data, rec["flat"], np.float32).reshape(-1, 2)
    unfolded = nib.load(unfold).agg_data("NIFTI_INTENT_POINTSET")
    np.testing.assert_allclose(flat[:, 0], -unfolded[:, 0])  # right hemi mirrored

    m = rec["metrics"]["thickness"]
    values = nib.load(metrics["thickness"]).agg_data()
    qv = decode(data, m["q"], np.uint8)
    np.testing.assert_allclose(
        m["qlo"] + qv / 254 * (m["qhi"] - m["qlo"]),
        values,
        atol=(values.max() - values.min()) / 254,
    )

    sf = rec["subfields"]
    assert {v[0] for v in sf["table"].values()} == {"CA1", "CA3"}
    assert set(decode(data, sf["q"], np.uint8)) == {1, 2}


def test_shared_limits_pool_subjects(surfaces):
    viewer = SurfaceViewer("Surfaces")
    for sub, native, unfold, metrics, labels in surfaces:
        viewer.add(sub, "L", "hipp", native, unfold, metrics)
    data = page_data(viewer.render())
    lo, hi = data["shared"]["hipp|thickness"]
    own = [r["metrics"]["thickness"]["own"] for r in data["records"]]
    assert lo < own[1][0] and hi > own[0][1]  # spans both subjects' ranges


def test_mismatched_surfaces_rejected(surfaces, tmp_path):
    sub, native, _, metrics, _ = surfaces[0]
    xy, faces = structured_sheet(11)
    other = tmp_path / "other.surf.gii"
    write_surface(other, np.c_[xy, np.zeros(len(xy))], faces)
    with pytest.raises(ValueError):
        SurfaceViewer("Surfaces").add(sub, "L", "hipp", native, other, metrics)


@pytest.mark.skipif(shutil.which("node") is None, reason="node not available")
def test_page_script_parses(surfaces, tmp_path):
    viewer = SurfaceViewer("Surfaces")
    for sub, native, unfold, metrics, labels in surfaces:
        viewer.add(sub, "L", "hipp", native, unfold, metrics, subfields=labels)
    page = viewer.render()
    script = page[page.rindex("<script>") + len("<script>") : page.rindex("</script>")]
    path = tmp_path / "viewer.js"
    path.write_text(script)
    subprocess.run(["node", "--check", str(path)], check=True)
