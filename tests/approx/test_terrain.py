

import numpy
import approx.terrain


def test_gen_knots():
    """
    This function tries to test generating of knot vector
    :return: None
    """
    knots = approx.terrain.gen_knots(7)
    assert numpy.array_equal(knots, numpy.array((0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0)))


# def test_base_factory():
#     """
#     This function tries to test basis_factory function
#     :return: None
#     """
#     knots = numpy.array((0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0))
#
#     num = 2
#     n_basf = len(knots) - 3
#     y_coords = [0.0] * n_basf
#     base_fact = approx.terrain.basis_factory(2)
#     for k in range(0, n_basf):
#         temp = [0.0] * (num + 1)
#         for i in range(0, num+1):
#             t_param = min(knots) + max(knots) * i / float(num)
#             temp[i] = base_fact(knots, k, t_param)
#         y_coords[k] = temp
#
#     assert y_coords == [[1.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]]


# def test_spline_base():
#     """
#     This function tries to test spline_base function
#     :return: None
#     """
#     knots = numpy.array((0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0))
#
#     num = 2
#     n_basf = len(knots) - 3
#     y_coords = [0.0] * n_basf
#     for k in range(0, n_basf):
#         temp = [0.0] * (num + 1)
#         for i in range(0, num+1):
#             t_param = min(knots) + max(knots) * i / float(num)
#             temp[i] = approx.terrain.spline_base(knots, k, t_param)
#         y_coords[k] = temp
#
#     assert y_coords == [[1.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]]
