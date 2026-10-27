from types import SimpleNamespace

import numpy as np
import pytest
from lib.qc import ERROR, WARNING, QCCheckError, QCResult, overall_status, verify
from lib.qc_report import QCReport
from lib.qc_srcsink import compute_qc, cyclic_runs, plot_loop, plot_surfaces
from lib.surface import assign_boundary_terminals


def make_result(passed_error=True, passed_warning=True):
    result = QCResult("demo", {"hemi": "L", "label": "hipp"})
    result.add_metric("n", np.int64(3), "a count")
    result.add_check("err", passed_error, "error check", severity=ERROR)
    result.add_check("warn", passed_warning, "warning check", severity=WARNING)
    return result


def test_status():
    assert make_result().status == "pass"
    assert make_result(passed_warning=False).status == "warn"
    assert make_result(passed_error=False).status == "fail"
    assert overall_status([make_result(), make_result(passed_warning=False)]) == "warn"


def test_invalid_severity():
    with pytest.raises(ValueError):
        QCResult("demo").add_check("x", True, severity="info")


def test_save_load_roundtrip(tmp_path):
    result = make_result(passed_warning=False)
    result.add_data("array", np.arange(3))
    path = tmp_path / "qc.json"
    result.save(path)

    loaded = QCResult.load(path)
    assert loaded.to_dict() == {**result.to_dict(), "data": {"array": [0, 1, 2]}}
    assert loaded.metrics["n"]["value"] == 3
    assert loaded.status == "warn"


def test_verify_raises_on_error_only():
    rows = verify([make_result(passed_warning=False)])
    assert len(rows) == 2

    with pytest.raises(QCCheckError, match="err"):
        verify([make_result(passed_error=False)], report="report.html")

    rows = verify([make_result(passed_error=False)], warn_only=True)
    assert [r["passed"] for r in rows] == [False, True]


def test_load_warn_checks(tmp_path):
    path = tmp_path / "qc.json"
    make_result(passed_error=False).save(path)

    assert QCResult.load(path).status == "fail"
    loaded = QCResult.load(path, warn_checks=["err"])
    assert loaded.status == "warn"
    rows = verify([loaded])
    assert [r["severity"] for r in rows] == [WARNING, WARNING]
    assert QCResult.load(path, warn_checks=["other"]).status == "fail"


def test_report_renders(tmp_path):
    results = [make_result(), make_result(passed_error=False)]
    report = QCReport("Demo", entities={"subject": "001"})
    report.add_summary(results)
    for r in results:
        section = report.add_section(r.label, status=r.status)
        section.add_checks(r)
        section.add_metrics(r)
    path = tmp_path / "report.html"
    report.save(path)

    page = path.read_text()
    assert page.startswith("<!doctype html>")
    assert "chip-fail" in page and "chip-pass" in page
    assert "prefers-color-scheme: dark" in page


# --- srcsink QC ---


def fan_disk(n):
    faces = np.array([[3, i, (i + 1) % n, n] for i in range(n)]).ravel()
    return SimpleNamespace(faces=faces)


def srcsink_qc(labels, hemi):
    n = len(labels)
    distances = 1.0 - np.eye(4)[labels]
    surface = fan_disk(n)
    out, info = assign_boundary_terminals(distances, np.arange(n), surface, 2)
    return compute_qc(distances, np.arange(n), out, info, hemi, {"hemi": hemi}, {})


def checks_by_name(result):
    return {c["name"]: c for c in result.checks}


def test_cyclic_runs():
    assert cyclic_runs([0, 0, 1, 1, 0]) == 2
    assert cyclic_runs([0, 0, 0]) == 1
    assert cyclic_runs([0, 2, 1, 3]) == 4


def test_srcsink_hemisphere_order():
    # fan disk winds with increasing index: this order is expected for hemi-L
    labels = np.repeat([0, 2, 1, 3], 6)
    assert checks_by_name(srcsink_qc(labels, "L"))["terminal_order_hemisphere"][
        "passed"
    ]
    check = checks_by_name(srcsink_qc(labels, "R"))["terminal_order_hemisphere"]
    assert not check["passed"]
    assert check["severity"] == ERROR


def test_srcsink_fragmentation_and_relabelling():
    labels = np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 4 + [2] * 3 + [3] * 5)
    checks = checks_by_name(srcsink_qc(labels, "L"))
    assert not checks["nearest_terminal_fragmentation"]["passed"]
    assert not checks["relabelled_fraction"]["passed"]  # 3 of 27 > 10%
    assert checks["terminal_order_hemisphere"]["passed"]
    assert all(c["severity"] == WARNING for c in checks.values() if not c["passed"])


def test_srcsink_missing_terminal():
    # AP sink absent from the nearest labels
    labels = np.repeat([0, 2, 3], 8)
    checks = checks_by_name(srcsink_qc(labels, "L"))
    assert not checks["nearest_terminal_fragmentation"]["passed"]
    assert "missing" in checks["nearest_terminal_fragmentation"]["message"]


def test_srcsink_figures():
    labels = np.repeat([0, 2, 1, 3], 6)
    result = srcsink_qc(labels, "L")
    n = len(labels)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    points = np.vstack(
        [np.c_[np.cos(theta), np.sin(theta), 0.1 * np.sin(2 * theta)], [0, 0, 0]]
    )
    faces = fan_disk(n).faces.reshape((-1, 4))[:, 1:]
    coords = np.linspace(0, 1, n + 1)

    assert plot_loop(result) is not None
    assert len(plot_surfaces(result, points, faces, coords, coords)) == 3
