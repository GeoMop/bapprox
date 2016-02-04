"""
This modules is used for testing public functions and methods from
module convert.bspline
"""

import convert.bspline


def test_approx_convertor_convert():
    """
    Test public method convert() form class BSplineSurfConvertor
    """
    scipy_bspline = [[0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                     [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                     [0.0, 0.1, 0.0, 0.1, 0.3, 0.1, 0.0, 0.1, 0.0],
                     3, 3]
    convertor = convert.bspline.BSplineSurfConvertor(scipy_bspline)
    occ_bspline = convertor.convert()
    assert occ_bspline is not None
