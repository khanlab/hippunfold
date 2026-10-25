import numpy as np
from lib.surface import isotonic_increasing, monotone_free_arcs


def test_isotonic_increasing():
    y = np.array([0.1, 0.3, 0.2, 0.4, 0.35, 0.36, 0.9])
    fit = isotonic_increasing(y)
    assert np.all(np.diff(fit) >= 0)
    assert np.isclose(fit.sum(), y.sum())  # pooling keeps the mean
    assert np.allclose(fit[[0, 6]], y[[0, 6]])
    assert np.allclose(isotonic_increasing([1, 2, 3]), [1, 2, 3])


def loop_setup():
    # loop of 12 vertices: src arc, free arc, sink arc, free arc (back to src)
    loop = np.arange(12)
    src = np.isin(np.arange(12), [0, 1, 2])
    sink = np.isin(np.arange(12), [6, 7, 8])
    values = np.zeros(12)
    values[sink] = 1
    values[[3, 4, 5]] = [0.3, 0.2, 0.7]  # not increasing from src to sink
    values[[9, 10, 11]] = [0.8, 0.6, 0.1]  # sink -> src, decreasing: fine
    return loop, values, src, sink


def test_monotone_free_arcs():
    loop, values, src, sink = loop_setup()
    fixed_values, n = monotone_free_arcs(loop, values, src, sink)
    assert n == 2  # vertices 3 and 4 pooled
    assert np.all(np.diff(fixed_values[[2, 3, 4, 5, 6]]) > 0)  # strictly
    assert np.all(np.diff(fixed_values[[8, 9, 10, 11]]) < 0)
    assert np.allclose(fixed_values[[9, 10, 11]], values[[9, 10, 11]], atol=1e-8)
    assert np.allclose(fixed_values[src | sink], values[src | sink])


def test_monotone_free_arcs_rotation_invariant():
    loop, values, src, sink = loop_setup()
    expected, _ = monotone_free_arcs(loop, values, src, sink)
    for shift in (1, 4, 7):
        rolled = np.roll(loop, shift)
        result, _ = monotone_free_arcs(rolled, values, src, sink)
        assert np.allclose(result, np.roll(expected, shift))
