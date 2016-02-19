"""
Simple test of terrain approximation
"""

import numpy
import approx
import approx.terrain
import convert
import convert.bspline

def test_terrain_approx():
    """
    This function is used for testing of terrain approximation
    """
    u_knots = approx.terrain.gen_knots(7)
    v_knots = approx.terrain.gen_knots(7)
    terrain_data = numpy.matrix([
        [0.0, 0.0, 0.0], [0.0, 0.5, 0.4], [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.2], [0.5, 0.5, 0.8], [0.5, 1.0, 0.2],
        [1.0, 0.0, 0.0], [1.0, 0.5, 0.4], [1.0, 1.0, 0.0]])
    return approx.terrain.approx(terrain_data, u_knots, v_knots)

if __name__ == '__main__':
    # test_spline_base()
    poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = test_terrain_approx()
    occ_bspline = convert.bspline.raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
