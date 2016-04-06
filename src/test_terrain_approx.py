"""
Simple test of terrain approximation
"""

import numpy
import time
import approx
import approx.terrain
import convert
import convert.bspline


def test_terrain_approx(method):
    """
    This function is used for testing of terrain approximation
    :param method: method used for approximation
    """
    u_knots = approx.terrain.gen_knots(7)
    v_knots = approx.terrain.gen_knots(7)
    # Show base functions for u knot vector
    # approx.terrain.test_spline_base(u_knots)
    terrain_data = numpy.matrix([
        [0.0, 0.0, 0.0, 1.0], [0.0, 0.5, 0.4, 1.0], [0.0, 1.0, 0.0, 1.0],
        [0.5, 0.0, 0.2, 1.0], [0.5, 0.5, 0.8, 1.0], [0.5, 1.0, 0.2, 1.0],
        [1.0, 0.0, 0.0, 1.0], [1.0, 0.5, 0.4, 1.0], [1.0, 1.0, 0.0, 1.0]])
    return approx.terrain.approx(method, terrain_data, u_knots, v_knots)


def display_results(occ_bspline, points):
    """
    Try to display OCC
    :param occ_bspline:
    :param points:
    :return:
    """

    # Initialize displaying
    from OCC.Display.SimpleGui import init_display
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.EraseAll()

    # Draw b-spline surfaces
    display.DisplayShape(occ_bspline.GetHandle(), update=True)

    # # Draw points of B-Spline
    import OCC.Prs3d
    import OCC.Quantity
    import OCC.Graphic3d
    import OCC.Aspect
    import OCC.gp

    a_presentation = OCC.Prs3d.Prs3d_Presentation(display._struc_mgr)
    group = OCC.Prs3d.Prs3d_Root_CurrentGroup(a_presentation.GetHandle()).GetObject()
    black = OCC.Quantity.Quantity_Color(OCC.Quantity.Quantity_NOC_BLACK)
    asp = OCC.Graphic3d.Graphic3d_AspectLine3d(black, OCC.Aspect.Aspect_TOL_SOLID, 1)

    gg = OCC.Graphic3d.Graphic3d_ArrayOfPoints(len(points),
                                               False,  # hasVColors
                                               )

    for point in points:
        pnt = OCC.gp.gp_Pnt(point[0], point[1], point[2])
        gg.AddVertex(pnt)

    group.SetPrimitivesAspect(asp.GetHandle())
    group.AddPrimitiveArray(gg.GetHandle())
    a_presentation.Display()

    start_display()


def test_svd_approx():
    """
    Test our approximation of B-Spline surface
    :return: None
    """
    # test_spline_base()
    poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = test_terrain_approx('svd')
    occ_bspline = convert.bspline.raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
    points = approx.terrain.gen_points(poles, u_knots, v_knots, u_mults, v_mults, u_num=50, v_num=50)
    display_results(occ_bspline, points)


def test_qr_approx():
    """
    Test our approximation of B-Spline surface
    :return: None
    """
    # test_spline_base()
    poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = test_terrain_approx('qr')
    occ_bspline = convert.bspline.raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
    points = approx.terrain.gen_points(poles, u_knots, v_knots, u_mults, v_mults, u_num=50, v_num=50)
    display_results(occ_bspline, points)


def test_scipy_approx():
    """
    Test SciPy approximation of B-Spline surface
    :return: None
    """
    terrain_data = [
        [0.0, 0.0, 0.0], [0.0, 0.5, 0.4], [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.2], [0.5, 0.5, 0.8], [0.5, 1.0, 0.2],
        [1.0, 0.0, 0.0], [1.0, 0.5, 0.4], [1.0, 1.0, 0.0]]
    tX = [item[0] for item in terrain_data]
    tY = [item[1] for item in terrain_data]
    tZ = [item[2] for item in terrain_data]
    from scipy import interpolate
    print('SciPy approximation ...')
    start_time = time.time()
    tck, fp, ior, msg = interpolate.bisplrep(tX, tY, tZ, kx=2, ky=2, full_output=1)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))
    occ_bspline = convert.bspline.scipy_to_occ(tck)
    # Compute difference between original terrain data and B-Spline surface
    u_num = v_num = 50
    points = [[float(i)/u_num, float(j)/u_num, 0.0] for i in range(v_num+1) for j in range(u_num+1)]
    points = [(it[0], it[1], interpolate.bisplev(it[0], it[1], tck)) for it in points]
    # points = terrain_data

    display_results(occ_bspline, points)


if __name__ == '__main__':
    test_svd_approx()
    # approx.terrain.test_spline_base_vec()
    # test_qr_approx()
    # test_scipy_approx()
    # u_knots = approx.terrain.gen_knots(9)
    # approx.terrain.test_spline_base_vec(u_knots)
    # approx.terrain.test_base_fact(u_knots)
    # approx.terrain.test_spline_base(u_knots)
    # # approx.terrain.norm_blend(u_knots, 1, 0.7)
    # approx.terrain.test_norm_blend(u_knots)
    # base_factory = approx.terrain.basis_factory(3)
    # t = 0.5
    # for i in range(3):
    #     out = base_factory(t, i, u_knots)
    #     print(i, out)
