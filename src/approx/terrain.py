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


def basis_factory(degree):
    """
    Returns a basis_function for the given degree
    :param degree:
    :return basis function
    """

    if degree == 0:

        def basis_function(knots, i, t):
            """
            The basis function for degree = 0
            :param knots: list of knot vectors
            :param i:
            :param t:
            """
            t_this = knots[i]
            t_next = knots[i+1]
            out = 1.0 if t_this <= t < t_next else 0.0
            return out
    else:

        def basis_function(knots, i, t):
            """
            The basis function for degree > 0
            :param knots: list of knots
            :param i:
            :param t:
            """
            out = 0.0
            try:
                t_this = knots[i]
                t_next = knots[i+1]
                t_precog = knots[i+degree]
                t_horizon = knots[i+degree+1]
            except IndexError:
                return 0.0

            top = t - t_this
            bottom = t_precog - t_this
            if bottom != 0:
                out = top / bottom * basis_factory(degree-1)(knots, i, t)

            top = t_horizon - t
            bottom = t_horizon - t_next
            if bottom != 0:
                out += top / bottom * basis_factory(degree-1)(knots, i+1, t)

            return out

    # "Enrich" the function with information about its "configuration"
    basis_function.lower = None if degree == 0 else basis_factory(degree-1)
    basis_function.degree = degree
    return basis_function


def test_base_fact(knots):
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
    base_fact = basis_factory(2)
    for k in range(0, n_basf):
        temp = {}
        for i in range(0, num+1):
            t_param = min(knots) + max(knots) * i / float(num)
            temp[i] = base_fact(knots, k, t_param)
            # temp[i] = norm_blend(knots, k, t_param)
        y_coords[k] = temp

    diff_x = (max(knots) - min(knots)) / num
    x_coord = [min(knots) + diff_x*i for i in range(0, num+1)]

    for temp in y_coords.values():
        plt.plot(x_coord, temp.values())
    # plt.show()


def spline_base(knot_vec, basis_fnc_idx, t_param):
    """
    This function compute normalized blending function aka base function of B-Spline curve or surface.
    This function implement some optimization.
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
        value = 0.0
        if knot_vec[basis_fnc_idx] <= t_param <= knot_vec[basis_fnc_idx+1] and kvn[basis_fnc_idx] != 0:
            value = (t_param - tks[0])**2 / ((tks[2] - tks[0]) * (tks[1] - tks[0]))
        elif knot_vec[basis_fnc_idx+1] <= t_param <= knot_vec[basis_fnc_idx+2] and kvn[basis_fnc_idx+1] != 0:
            value = ((t_param - tks[0]) * (tks[2] - t_param)) / ((tks[2] - tks[0]) * (tks[2] - tks[1])) + \
                   ((t_param - tks[1]) * (tks[3] - t_param)) / ((tks[3] - tks[1]) * (tks[2] - tks[1]))
        elif knot_vec[basis_fnc_idx+2] <= t_param <= knot_vec[basis_fnc_idx+3] and kvn[basis_fnc_idx+2] != 0:
            value = (tks[3] - t_param)**2 / ((tks[3] - tks[1]) * (tks[3] - tks[2]))
        SB_CACHE[(tuple(knot_vec), basis_fnc_idx, t_param)] = value

    return value


def spline_surface(poles, u_param, v_param, u_knots, v_knots, u_mults, v_mults):
    """
    Compute z coordinate of B-Surface u and v degree is 2
    :param poles: matrix of "poles"
    :param u_param: X coordinate in range <0, 1>
    :param v_param: Y coordinate in range <0, 1>
    :param u_knots: list of u knots
    :param v_knots: list of v knots
    :param u_mults: list of u multiplicities
    :param v_mults: list of v multiplicities
    :return: tuple of (x, y, z) coordinate of B-Spline surface
    """

    # "Decompress" knot vectors using multiplicities
    # e.g
    # u_knots: (0.0, 0.5, 1.0) u_mults: (3, 1, 3) will be converted to
    # _u_knot: (0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0)
    _u_knots = []
    _v_knots = []
    for idx,mult in enumerate(u_mults):
        _u_knots.extend([u_knots[idx]] * mult)
    for idx,mult in enumerate(v_mults):
        _v_knots.extend([v_knots[idx]] * mult)

    u_n_basf = len(_u_knots) - 3
    v_n_basf = len(_v_knots) - 3

    uf_mat = [0.0] * u_n_basf
    vf_mat = [0.0] * v_n_basf

    # Pre-compute base values of functions
    for k in range(0, u_n_basf):
        uf_mat[k] = spline_base(_u_knots, k, u_param)
    for k in range(0, v_n_basf):
        vf_mat[k] = spline_base(_v_knots, k, v_param)

    x_coord, y_coord, z_coord = 0.0, 0.0, 0.0

    # Compute point at B-Spline surface
    for i in range(0, u_n_basf):
        for j in range(0, v_n_basf):
            base_i_j = uf_mat[i] * vf_mat[j]
            x_coord += poles[i][j][0] * base_i_j
            y_coord += poles[i][j][1] * base_i_j
            z_coord += poles[i][j][2] * base_i_j

    return x_coord, y_coord, z_coord


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


def differences(terrain_data, poles, u_knots, v_knots, u_mults, v_mults):
    """
    Try to compute difference between terrain data and B-Spline surface approximation.
    :param terrain_data:
    :param poles:
    :param u_knots:
    :param v_knots:
    :return: list of differences
    """
    print('Computing differences ...')
    start_time = time.time()
    # Create list of values of differences
    diff = [0.0] * len(terrain_data)
    idx = 0
    for point in terrain_data:
        # X and Y coordinates should be equal to u and v parameters at
        # approximated points
        u_param = point[0, 0]
        v_param = point[0, 1]
        x_coord, y_coord, z_coord = spline_surface(poles, u_param, v_param, u_knots, v_knots, u_mults, v_mults)
        # assert x_coord == u_param
        # assert y_coord == v_param
        diff[idx] = abs(z_coord - point[0, 2])
        idx += 1
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))
    return diff


def gen_points(poles, u_knots, v_knots, u_mults, v_mults, u_num=100, v_num=100):
    """
    Generates points from B-Spline approximation
    :param poles:
    :param u_knots:
    :param v_knots:
    :param u_num:
    :param v_num:
    :return: list of points laying on B-Spline surface
    """
    # Create list of points on surface
    print('Computing points ...')
    start_time = time.time()
    points = [[0.0, 0.0, 0.0] for j in range(u_num * v_num)]
    for i in range(v_num):
        for j in range(u_num):
            idx = i * v_num + j
            u_param = float(i)/(u_num - 1)
            v_param = float(j)/(v_num - 1)
            x_coord, y_coord, z_coord = spline_surface(poles, u_param, v_param, u_knots, v_knots, u_mults, v_mults)
            points[idx] = (x_coord, y_coord, z_coord)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))
    return points


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

    g_z_mat = terrain_data[:, 2]

    print('Computing QR ...')
    start_time = time.time()
    q_mat, r_mat = numpy.linalg.qr(b_mat, mode='full')
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Computing Z matrix ...')
    start_time = time.time()
    z_mat = numpy.linalg.lstsq(r_mat, q_mat.transpose())[0] * g_z_mat
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

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

    return poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg
