import numpy as np
from lib.dentate_pd_labels import check_convention, swap_src_sink


def toy_dseg():
    """a dentate (8) strip inside grey matter (1), with the src (9) at the end
    exposed to the background and the sink (10) at the end buried in grey
    matter and SRLM (2)"""
    d = np.zeros((40, 20, 20), np.int16)
    d[5:35, 4:16, 4:16] = 1
    d[5:35, 8:12, 14:16] = 2
    d[8:30, 8:12, 8:12] = 8
    d[0:8, 8:12, 8:12] = 9  # reaches the edge of the grey matter
    d[30:33, 8:12, 8:12] = 10  # enclosed
    return d


def test_template_convention_detected():
    convention, details = check_convention(toy_dseg())
    assert convention == "template"
    assert details["src_neighbours"][0] > details["sink_neighbours"].get(0, 0)


def test_swapped_detected_and_fixed():
    swapped = swap_src_sink(toy_dseg())
    assert check_convention(swapped)[0] == "swapped"
    fixed = swap_src_sink(swapped)
    assert np.array_equal(fixed, toy_dseg())
    assert check_convention(fixed)[0] == "template"


def test_missing_labels():
    d = toy_dseg()
    d[d == 10] = 8
    assert check_convention(d)[0] == "missing"
