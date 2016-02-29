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
    u_knots = approx.terrain.gen_knots(9)
    v_knots = approx.terrain.gen_knots(9)
    # approx.terrain.test_spline_base(u_knots)
    terrain_data = numpy.matrix([
        [0.0, 0.0, 0.0], [0.0, 0.5, 0.4], [0.0, 1.0, 0.0],
        [0.5, 0.0, 0.2], [0.5, 0.5, 0.8], [0.5, 1.0, 0.2],
        [1.0, 0.0, 0.0], [1.0, 0.5, 0.4], [1.0, 1.0, 0.0]])
    return approx.terrain.approx(terrain_data, u_knots, v_knots)

if __name__ == '__main__':
    # test_spline_base()
    poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = test_terrain_approx()
    # poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg, points = test_terrain_approx()
    occ_bspline = convert.bspline.raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)

    # Initialize displaying
    from OCC.Display.SimpleGui import init_display
    # import OCC.Prs3d
    # import OCC.Quantity
    # import OCC.Graphic3d
    # import OCC.Aspect
    # import OCC.gp
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.EraseAll()

    # Draw b-spline surfaces
    display.DisplayShape(occ_bspline.GetHandle(), update=True)

    # Draw points of B-Spline
    # a_presentation = OCC.Prs3d.Prs3d_Presentation(display._struc_mgr)
    # group = OCC.Prs3d.Prs3d_Root_CurrentGroup(a_presentation.GetHandle()).GetObject()
    # black = OCC.Quantity.Quantity_Color(OCC.Quantity.Quantity_NOC_BLACK)
    # asp = OCC.Graphic3d.Graphic3d_AspectLine3d(black, OCC.Aspect.Aspect_TOL_SOLID, 1)
    #
    # gg = OCC.Graphic3d.Graphic3d_ArrayOfPoints(len(points),
    #                                            False,  # hasVColors
    #                                            )
    #
    # for point in points:
    #     pnt = OCC.gp.gp_Pnt(point[0], point[1], point[2])
    #     gg.AddVertex(pnt)
    #
    # group.SetPrimitivesAspect(asp.GetHandle())
    # group.AddPrimitiveArray(gg.GetHandle())
    # a_presentation.Display()

    start_display()
