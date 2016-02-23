"""
This module tries to approximate 2.5D array of terrain points
using B-Spline surface.
"""

import matplotlib.pyplot as plt
import numpy
import numpy.linalg

__author__ = 'Jiri Hnidek <jiri.hnidek@tul.cz>, Jiri Kopal <jiri.kopal@tul.cz>'

# Cache used of knot vector (computation of differences)
KVN_CACHE = {}


def spline_base(knot_vec, basis_fnc_idx, t_param):
    """
    :param knot_vec: knot vector
    :param basis_fnc_idx: index of basis function
    :param t_param: parameter t in interval <0, 1>
    :return: value of basis function
    """

    # When basis function has zero value at given interval, then return 0
    if t_param < knot_vec[basis_fnc_idx] or knot_vec[basis_fnc_idx+3] < t_param:
        return 0.0

    try:
        kvn = KVN_CACHE[tuple(knot_vec)]
    except KeyError:
        knot_vec_len = len(knot_vec)
        kvn = [0] * knot_vec_len
        i = 0
        while i < knot_vec_len - 1:
            if knot_vec[i] - knot_vec[i+1] != 0:
                kvn[i] = 1.0
            i += 1
        KVN_CACHE[tuple(knot_vec)] = kvn

    tks = [knot_vec[basis_fnc_idx + k] for k in range(0, 4)]

    if knot_vec[basis_fnc_idx] <= t_param <= knot_vec[basis_fnc_idx+1] and kvn[basis_fnc_idx] != 0:
        return (t_param - tks[0])**2 / ((tks[2] - tks[0]) * (tks[1] - tks[0]))
    elif knot_vec[basis_fnc_idx+1] <= t_param <= knot_vec[basis_fnc_idx+2] and kvn[basis_fnc_idx+1] != 0:
        return ((t_param - tks[0]) * (tks[2] - t_param)) / ((tks[2] - tks[0]) * (tks[2] - tks[1])) + \
               ((t_param - tks[1]) * (tks[3] - t_param)) / ((tks[3] - tks[1]) * (tks[2] - tks[1]))
    elif knot_vec[basis_fnc_idx+2] <= t_param <= knot_vec[basis_fnc_idx+3] and kvn[basis_fnc_idx+2] != 0:
        return (tks[3] - t_param)**2 / ((tks[3] - tks[1]) * (tks[3] - tks[2]))

    assert False, "Should not reach this point"


def test_spline_base():
    """
    This function tries to test  spline_base function
    :return: None
    """
    # knots = [0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]
    # knots = numpy.array((0.0, 0.0, 0.0, 1/3.0, 2/3.0, 1.0, 1.0, 1.0))
    knots = gen_knots(7)
    print(knots)
    num = 100
    n_basf = len(knots) - 3
    y_coords = {}
    for k in range(0, n_basf):
        temp = {}
        for i in range(0, num+1):
            t_param = min(knots) + max(knots) * i / float(num)
            temp[i] = spline_base(knots, k, t_param)
        y_coords[k] = temp

    diff_x = (max(knots) - min(knots)) / num
    x_coord = [min(knots) + diff_x*i for i in range(0, num+1)]

    for temp in y_coords.values():
        plt.plot(x_coord, temp.values())
    plt.show()


def gen_knots(num=10):
    """
    This function generates vector of knots according number
    :param num: length of vector
    :return: array of knots
    """
    knots = numpy.array((0.0,) * num)
    diff = 1.0 / (num - 3 - 3 + 1)
    for i in range(3, num - 3 + 1):
        knots[i] = (i - 2) * diff
    knots[-3:] = 1.0
    return knots


def transform(terrain_data):
    """
    This function will transform x and y coordinates to range <0.0, 1.0>
    :param terrain_data: matrix of 3D terrain data
    :return: transformed matrix of 3D terrain data
    """



def approx(terrain_data, u_knots, v_knots):
    """
    This function tries to approximate terrain data with B-Spline surface patches
    :param terrain_data: matrix of 3D terrain data
    :param u_knots: array of u knots
    :param v_knots: array of v knots
    :return: B-Spline patch
    """
    num_pnt = terrain_data.shape[0]

    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3

    b_mat = numpy.zeros((num_pnt, u_n_basf * v_n_basf))

    uf_mat = numpy.matrix((0.0,) * u_n_basf)
    vf_mat = numpy.matrix((0.0,) * v_n_basf)

    # Black magic alias linear algebra :-)
    for j in range(0, num_pnt):
        for k in range(0, u_n_basf):
            uf_mat[(0, k)] = spline_base(u_knots, k, terrain_data[j, 0])
        for k in range(0, v_n_basf):
            vf_mat[(0, k)] = spline_base(v_knots, k, terrain_data[j, 1])
        b_mat[j] = vf_mat * numpy.kron(numpy.eye(v_n_basf), uf_mat)

    g_mat = terrain_data[:, 2]

    q_mat, r_mat = numpy.linalg.qr(b_mat, mode='full')

    z_mat = numpy.linalg.lstsq(r_mat, q_mat.transpose())[0] * g_mat

    # Create list of poles from z_mat
    poles = [[[0.0, 0.0, 0.0] for i in range(v_n_basf)] for j in range(u_n_basf)]
    for i in range(0, u_n_basf):
        for j in range(0, v_n_basf):
            x_coord = float(i)/(u_n_basf - 1)
            y_coord = float(j)/(v_n_basf - 1)
            z_coord = z_mat[i * v_n_basf + j, 0]
            poles[i][j] = (x_coord, y_coord, z_coord)

    # Create degrees
    u_deg = 2
    v_deg = 2

    # Convert knot vectors
    u_knots = list(set(u_knots))
    u_knots.sort()
    v_knots = list(set(v_knots))
    v_knots.sort()

    # Create vectors of multiplicities
    u_mults = [1] * (u_n_basf - 1)
    u_mults[0] = u_mults[-1] = 3
    v_mults = [1] * (v_n_basf - 1)
    v_mults[0] = v_mults[-1] = 3

    return poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg
