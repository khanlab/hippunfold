from types import SimpleNamespace

import networkx as nx
import numpy as np
import pytest
from lib.surface import (
    TERMINAL_ORDERS,
    assign_boundary_terminals as assign_with_info,
    get_boundary_loop,
)


def assign_boundary_terminals(*args, **kwargs):
    labels, _ = assign_with_info(*args, **kwargs)
    return labels


def fan_disk(n):
    """Triangulated disk: boundary loop 0..n-1, centre vertex n."""
    faces = np.array([[3, i, (i + 1) % n, n] for i in range(n)]).ravel()
    return SimpleNamespace(faces=faces)


def annulus_disk(n):
    """
    Triangulated disk with a band: outer loop 0..n-1, inner ring n..2n-1,
    centre vertex 2n. Outer loop and inner ring are the boundary band.
    """
    faces = []
    for i in range(n):
        j = (i + 1) % n
        faces += [[3, i, j, n + i], [3, j, n + j, n + i], [3, n + i, n + j, 2 * n]]
    return SimpleNamespace(faces=np.array(faces).ravel())


def one_hot_distances(labels):
    """Distance 0 to a vertex's own terminal and 1 to the others."""
    return 1.0 - np.eye(4)[labels]


def cyclic_runs(labels):
    """Labels of the runs around a cyclic sequence."""
    runs = [lab for i, lab in enumerate(labels) if i == 0 or lab != labels[i - 1]]
    if len(runs) > 1 and runs[0] == runs[-1]:
        runs.pop()
    return runs


def assert_valid_order(labels):
    runs = cyclic_runs(list(labels))
    assert len(runs) == 4, runs
    start = runs.index(0)
    assert tuple(runs[start:] + runs[:start]) in TERMINAL_ORDERS, runs


def assign(labels, surface=None, nmin=2, **kwargs):
    n = len(labels)
    surface = fan_disk(n) if surface is None else surface
    return assign_boundary_terminals(
        one_hot_distances(np.asarray(labels)), np.arange(n), surface, nmin, **kwargs
    )


def test_valid_labels_unchanged():
    labels = np.repeat([0, 2, 1, 3], 6)
    np.testing.assert_array_equal(assign(labels), labels)


def test_mirrored_order_unchanged():
    labels = np.repeat([0, 3, 1, 2], 6)
    np.testing.assert_array_equal(assign(labels), labels)


def test_multi_vertex_stray_is_merged():
    # example from issue #580: PD src (2) has a 3-vertex stray inside PD sink (3)
    labels = np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 4 + [2] * 3 + [3] * 5)
    expected = np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 12)
    np.testing.assert_array_equal(assign(labels, nmin=5), expected)


def test_adjacent_ap_terminals_are_separated():
    # AP src and sink next to each other, so PD terminals are not opposite
    labels = np.repeat([0, 1, 2, 3], 6)
    result = assign(labels)
    assert_valid_order(result)
    # cheapest fix: insert an nmin-long PD arc between the AP terminals, and
    # absorb one of the (6-vertex) PD arcs on the other side
    assert np.sum(result != labels) == 2 + 6


def test_order_violation_via_band():
    # PD sink appears on both sides of AP src (as seen in ds002168 L dentate);
    # each label is contiguous on the band, but not on the loop
    n = 30
    loop_labels = np.array([2] * 10 + [3] * 3 + [0] * 3 + [3] * 10 + [1] * 4)
    labels = np.concatenate([loop_labels, loop_labels])  # inner ring matches
    surface = annulus_disk(n)

    result = assign_boundary_terminals(
        one_hot_distances(labels), np.arange(2 * n), surface, 2
    )

    assert_valid_order(result[:n])


def test_band_labels_are_contiguous():
    n = 30
    loop_labels = np.repeat([0, 2, 1, 3], [6, 9, 6, 9])
    rng = np.random.default_rng(0)
    # noisy inner ring, which should follow the loop
    labels = np.concatenate([loop_labels, rng.integers(0, 4, n)])
    surface = annulus_disk(n)

    result = assign_boundary_terminals(
        one_hot_distances(labels), np.arange(2 * n), surface, 2
    )

    np.testing.assert_array_equal(result[:n], loop_labels)
    faces = surface.faces.reshape((-1, 4))[:, 1:]
    band = nx.Graph(
        (u, v)
        for f in faces
        for u, v in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0]))
        if u < 2 * n and v < 2 * n
    )
    for k in range(4):
        nodes = np.flatnonzero(result == k)
        assert nx.number_connected_components(band.subgraph(nodes)) == 1


def test_min_terminal_vertices():
    # AP sink has only one vertex
    labels = np.array([0] * 6 + [2] * 9 + [1] * 1 + [3] * 9)
    result = assign(labels, nmin=4)
    assert_valid_order(result)
    assert np.bincount(result, minlength=4).min() >= 4


def test_loop_too_short():
    with pytest.raises(ValueError, match="too few"):
        assign(np.repeat([0, 2, 1, 3], 2), nmin=3)


def test_ap_length_penalty_shortens_ap():
    n = 40
    labels = np.repeat([0, 2, 1, 3], 10)
    # graded distances: AP src/sink nearly as close as PD over their arcs
    distances = one_hot_distances(labels) * 0.1
    surface = fan_disk(n)

    plain = assign_boundary_terminals(distances, np.arange(n), surface, 3)
    penalised = assign_boundary_terminals(
        distances, np.arange(n), surface, 3, ap_length_penalty=1.0
    )

    np.testing.assert_array_equal(plain, labels)
    assert_valid_order(penalised)
    assert np.sum(penalised == 0) == 3
    assert np.sum(penalised == 1) == 3


def test_boundary_vertex_ids_need_not_be_contiguous():
    labels = np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 4 + [2] * 3 + [3] * 5)
    n = len(labels)
    boundary_vertices = np.arange(n) * 3 + 7
    remap = dict(zip(range(n), boundary_vertices))
    remap[n] = 1000
    faces = fan_disk(n).faces.reshape((-1, 4))
    faces[:, 1:] = np.vectorize(remap.get)(faces[:, 1:])

    result = assign_boundary_terminals(
        one_hot_distances(labels),
        boundary_vertices,
        SimpleNamespace(faces=faces.ravel()),
        5,
    )

    np.testing.assert_array_equal(
        result, np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 12)
    )


def test_boundary_loop_follows_face_winding():
    n = 12
    faces = fan_disk(n).faces.reshape((-1, 4))[:, 1:]
    loop = get_boundary_loop(faces, np.arange(n))
    start = int(np.flatnonzero(loop == 0)[0])
    np.testing.assert_array_equal(np.roll(loop, -start), np.arange(n))


def test_info_reports_order_and_nearest_labels():
    labels = np.array([0] * 5 + [2] * 5 + [1] * 5 + [3] * 4 + [2] * 3 + [3] * 5)
    n = len(labels)

    result, info = assign_with_info(
        one_hot_distances(labels), np.arange(n), fan_disk(n), 5
    )

    # fan disk faces wind with increasing vertex index
    loop = info["loop"]
    start = int(np.flatnonzero(loop == 0)[0])
    np.testing.assert_array_equal(np.roll(loop, -start), np.arange(n))
    assert info["order"] == TERMINAL_ORDERS[0]
    np.testing.assert_array_equal(info["loop_nearest_labels"], labels[loop])
    np.testing.assert_array_equal(info["loop_labels"], result[loop])
    assert info["n_unreached"] == 0

    # reversing the winding mirrors the order
    faces = fan_disk(n).faces.reshape((-1, 4))
    faces[:, [1, 2]] = faces[:, [2, 1]]
    _, info = assign_with_info(
        one_hot_distances(labels), np.arange(n), SimpleNamespace(faces=faces.ravel()), 5
    )
    assert info["order"] == TERMINAL_ORDERS[1]
