"""
This module tries to approximate 2.5D array of terrain points
using B-Spline surface.
"""

import matplotlib.pyplot as plt
import numpy
import time
import numpy.linalg
# import scipy.sparse

__author__ = 'Jiri Hnidek <jiri.hnidek@tul.cz>, Jiri Kopal <jiri.kopal@tul.cz>'

# Cache used of knot vector (computation of differences)
KVN_CACHE = {}
SB_CACHE = {}


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
        value = SB_CACHE[(tuple(knot_vec), basis_fnc_idx, t_param)]
    except KeyError:
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
            value = (t_param - tks[0])**2 / ((tks[2] - tks[0]) * (tks[1] - tks[0]))
        elif knot_vec[basis_fnc_idx+1] <= t_param <= knot_vec[basis_fnc_idx+2] and kvn[basis_fnc_idx+1] != 0:
            value = ((t_param - tks[0]) * (tks[2] - t_param)) / ((tks[2] - tks[0]) * (tks[2] - tks[1])) + \
                   ((t_param - tks[1]) * (tks[3] - t_param)) / ((tks[3] - tks[1]) * (tks[2] - tks[1]))
        elif knot_vec[basis_fnc_idx+2] <= t_param <= knot_vec[basis_fnc_idx+3] and kvn[basis_fnc_idx+2] != 0:
            value = (tks[3] - t_param)**2 / ((tks[3] - tks[1]) * (tks[3] - tks[2]))
        SB_CACHE[(tuple(knot_vec), basis_fnc_idx, t_param)] = value

    return value


def spline_surface(x_coord, y_coord, u_knots, v_knots, z_mat):
    """
    Compute z coordinate of surface
    :param x_coord: X coordinate in range <0, 1>
    :param y_coord: Y coordinate in range <0, 1>
    :param u_knots: list of u knots
    :param v_knots: list of v knots
    :param z_mat: matrix of "poles"
    :return: Z coordinate of B-Spline surface
    """
    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3

    uf_mat = numpy.matrix((0.0,) * u_n_basf)
    vf_mat = numpy.matrix((0.0,) * v_n_basf)

    eye_mat = numpy.eye(v_n_basf)
    for k in range(0, u_n_basf):
        uf_mat[(0, k)] = spline_base(u_knots, k, x_coord)
    for k in range(0, v_n_basf):
        vf_mat[(0, k)] = spline_base(v_knots, k, y_coord)
    # print('vf_mat:', vf_mat)
    # print('uv_mat:', uf_mat)
    # print('eye:', eye_mat)
    # print('z_mat:', z_mat)
    z_coord = vf_mat * numpy.kron(eye_mat, uf_mat) * z_mat

    return z_coord[0, 0]


def test_spline_base(knots):
    """
    This function tries to test  spline_base function
    :param knots NumPy array of knots
    :return: None
    """
    # knots = numpy.array((0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0))
    # knots = numpy.array((0.0, 0.0, 0.0, 1/3.0, 2/3.0, 1.0, 1.0, 1.0))
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

    # Own computation of approximation
    print('Creating B matrix ...')
    start_time = time.time()
    eye_mat = numpy.eye(v_n_basf)
    # eye_mat = scipy.sparse.eye(v_n_basf)
    for j in range(0, num_pnt):
        for k in range(0, u_n_basf):
            uf_mat[(0, k)] = spline_base(u_knots, k, terrain_data[j, 0])
        for k in range(0, v_n_basf):
            vf_mat[(0, k)] = spline_base(v_knots, k, terrain_data[j, 1])
        b_mat[j] = vf_mat * numpy.kron(eye_mat, uf_mat)
        # b_mat[j] = vf_mat * scipy.sparse.kron(eye_mat, uf_mat)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    g_mat = terrain_data[:, 2]

    print('Computing QR ...')
    start_time = time.time()
    q_mat, r_mat = numpy.linalg.qr(b_mat, mode='full')
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Computing Z matrix ...')
    start_time = time.time()
    z_mat = numpy.linalg.lstsq(r_mat, q_mat.transpose())[0] * g_mat
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    # Create list of differences between terrain and points on surface
    print('Computing differences ...')
    start_time = time.time()
    tW = [0.0 for it in terrain_data]
    idx = 0
    for point in terrain_data:
        x_coord = point[0, 0]
        y_coord = point[0, 1]
        z_coord = spline_surface(x_coord, y_coord, u_knots, v_knots, z_mat)
        tW[idx] = abs(z_coord - point[0, 2])
        idx += 1
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    # Create list of points on surface
    # u_num = v_num = 2
    # points = [[float(i)/u_num, float(j)/u_num, 0.0] for i in range(v_num+1) for j in range(u_num+1)]
    # for i in range(v_num+1):
    #     for j in range(u_num+1):
    #         idx = i * v_num + j
    #         x_coord = points[idx][0]
    #         y_coord = points[idx][1]
    #         z_coord = spline_surface(x_coord, y_coord, u_knots, v_knots, z_mat)
    #         points[idx][2] = z_coord

    # Create list of poles from z_mat
    poles = [[[0.0, 0.0, 0.0] for i in range(v_n_basf)] for j in range(u_n_basf)]
    for i in range(0, u_n_basf):
        for j in range(0, v_n_basf):
            x_coord = float(i)/(u_n_basf - 1)
            y_coord = float(j)/(v_n_basf - 1)
            z_coord = z_mat[i * v_n_basf + j, 0]
            # For some reason we have to switch x and y coordinates.
            # Shame on me, but I don't know why :-/
            poles[i][j] = (y_coord, x_coord, z_coord)

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

    return poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg, tW
