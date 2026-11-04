import numpy as np
from lib.myelin import divide_t1_by_t2


def test_numerator_is_t1():
    np.testing.assert_allclose(divide_t1_by_t2([4.0, 9.0], [2.0, 3.0]), [2.0, 3.0])


def test_divide_by_zero():
    out = divide_t1_by_t2([4.0, -4.0, 0.0], [0.0, 0.0, 0.0])
    np.testing.assert_array_equal(out, [1000.0, -1000.0, 0.0])
    assert np.all(np.isfinite(out))


def test_integer_inputs_are_not_truncated():
    out = divide_t1_by_t2(np.array([3], dtype=np.int16), np.array([2], dtype=np.int16))
    np.testing.assert_allclose(out, [1.5])
