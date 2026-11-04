import numpy as np


def divide_t1_by_t2(t1, t2, zero_div_value=1000.0):
    """Voxelwise T1w/T2w ratio (myelin proxy).

    Division by zero gives +/-zero_div_value (sign of the numerator), and 0/0 gives 0,
    matching the intent of the previous `c3d ... -divide -replace inf 1000 -inf -1000
    NaN 0` command.

    Done in numpy rather than with c3d -divide, since c3d's argument order is the
    reverse of its documentation (`c3d A B -divide` gives B/A, not A/B), and its
    divide-by-zero results (DBL_MAX, not inf/NaN) are not caught by -replace and end
    up as +inf in the float32 output.
    """
    t1 = np.asarray(t1, dtype=np.float64)
    t2 = np.asarray(t2, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = t1 / t2
    return np.nan_to_num(ratio, nan=0.0, posinf=zero_div_value, neginf=-zero_div_value)
