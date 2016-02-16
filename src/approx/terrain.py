"""
This module tries to approximate 2.5D array of terrain points
using B-Spline surface.
"""

__author__ = 'Jiri Hnidek <jiri.hnidek@tul.cz>, Jiri Kopal <jiri.kopal@tul.cz>'


import matplotlib.pyplot as plt


def spline_base(knot_vec, basis_fnc_idx, t_param):
    """
    :param knot_vec:
    :param basis_fnc_idx:
    :param t_param:
    :return:
    """

    # When basis function has zero value at given interval, then return 0
    if t_param < knot_vec[basis_fnc_idx] or knot_vec[basis_fnc_idx+3] < t_param:
        return 0.0

    knot_vec_len = len(knot_vec)
    N = [0] * knot_vec_len

    i = 0
    while i < knot_vec_len - 1:
        if knot_vec[i] - knot_vec[i+1] != 0:
            N[i] = 1.0
        i += 1

    tks = [knot_vec[basis_fnc_idx + k] for k in range(0, 4)]

    if knot_vec[basis_fnc_idx] <= t_param <= knot_vec[basis_fnc_idx+1] and N[basis_fnc_idx] != 0:
        return (t_param - tks[0])**2 / ((tks[2] - tks[0]) * (tks[1] - tks[0]))
    elif knot_vec[basis_fnc_idx+1] <= t_param <= knot_vec[basis_fnc_idx+2] and N[basis_fnc_idx+1] != 0:
        return ((t_param - tks[0]) * (tks[2] - t_param)) / ((tks[2] - tks[0]) * (tks[2] - tks[1])) + \
               ((t_param - tks[1]) * (tks[3] - t_param)) / ((tks[3] - tks[1]) * (tks[2] - tks[1]))
    elif knot_vec[basis_fnc_idx+2] <= t_param <= knot_vec[basis_fnc_idx+3] and N[basis_fnc_idx+2] != 0:
        return (tks[3] - t_param)**2 / ((tks[3] - tks[1]) * (tks[3] - tks[2]))

    assert False, "Should not reach this point"


def test_spline_base():
    """
    This function tries to test  spline_base function
    :return: None
    """
    # knots = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    knots = [0.0, 0.0, 0.0, 1/3.0, 2/3.0, 1.0, 1.0, 1.0]
    num = 100
    n_basf = len(knots) - 3
    y = {}
    for k in range(0, n_basf):
        temp = {}
        for i in range(0, num+1):
            t_param = min(knots) + max(knots) * i / float(num)
            temp[i] = spline_base(knots, k, t_param)
        y[k] = temp

    dx = (max(knots) - min(knots)) / num
    x = [min(knots) + dx*i for i in range(0, num+1)]

    for temp in y.values():
        plt.plot(x, temp.values())
    plt.show()

if __name__ == '__main__':
    test_spline_base()