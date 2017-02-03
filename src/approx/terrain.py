"""
This module tries to approximate 2.5D array of terrain points
using B-Spline surface.
"""

import matplotlib.pyplot as plt
import numpy
import math
import time
import numpy.linalg
import scipy.sparse
import scipy.sparse.linalg

__author__ = 'Jiri Hnidek <jiri.hnidek@tul.cz>, Jiri Kopal <jiri.kopal@tul.cz>'

# Cache used of knot vector (computation of differences)
KVN_CACHE = {}
SB_CACHE = {}


# TODO: this function contains small bug. It return 0.0 instead 1.0 for t_param = 1.0
# This function is not used for production, but it is here only for testing.
def basis_factory(degree):
    """
    Returns a basis_function for the given degree
    :param degree:
    :return basis function
    """

    if degree == 0:
        def basis_function(knots, idx, t_param):
            """
            The basis function for degree = 0
            :param knots: list of knot vectors
            :param idx:
            :param t_param:
            """
            t_this = knots[idx]
            t_next = knots[idx + 1]
            out = 1.0 if t_this <= t_param < t_next else 0.0
            return out
    else:
        def basis_function(knots, idx, t_param):
            """
            The basis function for degree > 0
            :param knots: list of knots
            :param idx:
            :param t_param:
            """
            out = 0.0
            try:
                t_this = knots[idx]
                t_precog = knots[idx + degree]
            except IndexError:
                pass
            else:
                top = t_param - t_this
                bottom = t_precog - t_this
                if bottom != 0:
                    out = top / bottom * basis_factory(degree - 1)(knots, idx, t_param)
            try:
                t_next = knots[idx + 1]
                t_horizon = knots[idx + degree + 1]
            except IndexError:
                pass
            else:
                top = t_horizon - t_param
                bottom = t_horizon - t_next
                if bottom != 0:
                    out += top / bottom * basis_factory(degree - 1)(knots, idx + 1, t_param)

            return out

    # "Enrich" the function with information about its "configuration"
    basis_function.lower = None if degree == 0 else basis_factory(degree - 1)
    basis_function.degree = degree
    return basis_function


def spline_base_vec(knot_vec, t_param, order, sparse=False):
    """
    This function compute normalized blending function aka base function of B-Spline curve or surface.
    :param knot_vec:
    :param t_param:
    :param order: (0: function value, 1: derivative function value)
    :param sparse:
    :return:
    """

    def find_index(_knot_vec, _t_param):
        """
        This function try to find index for given t_param in knot_vec that
        is covered by all (3) base functions.
        :param _knot_vec:
        :param _t_param:
        :return:
        """

        n = len(_knot_vec)
        mn = 2
        mx = n - 3
        est = min((2 + int(_t_param * (n - 5))), n - 4)
        _idx = 0

        if _t_param >= _knot_vec[est]:
            mn = est
        elif _t_param < _knot_vec[est]:
            mx = est

        s = int(max(math.ceil(math.log(mx - mn, 2)), 1))

        for p in range(s + 1):
            if _t_param < _knot_vec[mn + 1]:
                _idx = mn - 2
                break
            elif _t_param > _knot_vec[mx - 1]:
                _idx = mx - 3
                break

            mid = mn + int((mx - mn) / 2)
            if mid != mn:
                if _t_param <= _knot_vec[mid]:
                    mx = mid
                elif _t_param > _knot_vec[mid]:
                    mn = mid
            else:
                _idx = mn - 2
                break

        return _idx

    idx = find_index(knot_vec, t_param)
    n_basf = len(knot_vec) - 3

    # Create sparse matrix
    if sparse is True:
        basis_values = numpy.zeros(3)
    else:
        basis_values = numpy.zeros(n_basf)

    tk1 = knot_vec[idx + 1]
    tk2 = knot_vec[idx + 2]
    tk3 = knot_vec[idx + 3]
    tk4 = knot_vec[idx + 4]

    d31 = tk3 - tk1
    d32 = tk3 - tk2
    d42 = tk4 - tk2

    dt1 = t_param - tk1
    dt2 = t_param - tk2
    d3t = tk3 - t_param
    d4t = tk4 - t_param

    d31_d32 = d31 * d32
    d42_d32 = d42 * d32

    # basis function values
    if order == 0:
        if sparse is True:
            basis_values[0] = (d3t * d3t) / d31_d32
            basis_values[1] = ((dt1 * d3t) / d31_d32) + ((dt2 * d4t) / d42_d32)
            basis_values[2] = (dt2 * dt2) / d42_d32

        else:
            basis_values[idx] = (d3t * d3t) / d31_d32
            basis_values[idx + 1] = ((dt1 * d3t) / d31_d32) + ((dt2 * d4t) / d42_d32)
            basis_values[idx + 2] = (dt2 * dt2) / d42_d32

    # basis function derivatives
    elif order == 1:
        if sparse is True:
            basis_values[0] = -2 * d3t / d31_d32
            basis_values[1] = (d3t - dt1) / d31_d32 + (d4t - dt2) / d42_d32
            basis_values[2] = 2 * dt2 / d42_d32
        else:
            basis_values[idx] = -2*d3t / d31_d32
            basis_values[idx + 1] = (d3t - dt1) / d31_d32 + (d4t - dt2) / d42_d32
            basis_values[idx + 2] = 2 * dt2 / d42_d32

    return basis_values, idx


def test_spline_base_vec(knots=numpy.array((0.0, 0.0, 0.0, 1/3.0, 2/3.0, 1.0, 1.0, 1.0)), sparse=False):
    """
    Test optimized version of spline base function
    :param knots: numpy array of knots
    :param sparse: is sparse matrix used
    :return:
    """

    num = 100
    n_basf = len(knots) - 3
    y_coords = {}
    for k in range(0, n_basf):
        temp = {}
        for i in range(0, num+1):
            t_param = min(knots) + max(knots) * i / float(num)
            if sparse is True:
                temp[i] = spline_base_vec(knots, t_param, 0, sparse)[0].toarray()[0]
            else:
                temp[i] = spline_base_vec(knots, t_param, 0, sparse)[0]
        y_coords[k] = temp

    diff_x = (max(knots) - min(knots)) / num
    x_coord = [min(knots) + diff_x*i for i in range(0, num+1)]

    for temp in y_coords.values():
        plt.plot(x_coord, temp.values())
    plt.show()


def build_ls_matrix(u_knots, v_knots, terrain, sparse):
    """
    Construction of the matrix (B) of the system of linear algebraic
    equations for control points of the 2th order B-spline surface
    :param u_knots:
    :param v_knots:
    :param terrain:
    :param sparse:
    :return:
    """
    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3
    terrain_len = len(terrain)

    row = numpy.zeros(terrain_len * 9)
    col = numpy.zeros(terrain_len * 9)
    data = numpy.zeros(terrain_len * 9)

    if sparse is True:
        nnz_b = 0
    else:
        mat_b = numpy.empty((terrain_len, u_n_basf * v_n_basf))

    interval = numpy.empty((terrain_len, 2))

    for idx in range(0, terrain_len):
        u_base_vec, i_idx = spline_base_vec(u_knots, terrain[idx, 0], 0, sparse)
        v_base_vec, j_idx = spline_base_vec(v_knots, terrain[idx, 1], 0, sparse)
        if sparse is True:
            # Hard-coded Kronecker product (problem based)
            for n in range(0, 3):
                data[nnz_b + 3 * n:nnz_b + 3 * (n + 1)] = v_base_vec[n] * u_base_vec
                for m in range(0, 3):
                    col[nnz_b + (3 * n) + m] = (j_idx + n) * u_n_basf + i_idx + m
            row[nnz_b:nnz_b+9] = idx
            nnz_b += 9
        else:
            mat_b[idx] = numpy.kron(v_base_vec, u_base_vec)

        interval[idx][0] = i_idx
        interval[idx][1] = j_idx

    if sparse is True:
        mat_b = scipy.sparse.csr_matrix((data, (row, col)), shape=(terrain_len, u_n_basf * v_n_basf))

    return mat_b, interval


def build_sparse_reg_matrix(u_knots, v_knots, quad):
    """
    Construction of the regularization matrix (A) to decrease variation of the terrain
    B z = b  ---> (B^T B + A)z = B^T b
    :param u_knots: vector of v-knots
    :param v_knots: vector of u-knots
    :param quad: points defining quadrangle area (array)
    :return: matrix
    """

    a = quad[:, 3] - quad[:, 2]
    b = quad[:, 0] - quad[:, 1]
    c = quad[:, 1] - quad[:, 2]
    d = quad[:, 0] - quad[:, 3]

    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3
    u_n_inter = len(u_knots) - 5
    v_n_inter = len(v_knots) - 5

    q_points = [0, (0.5 - 1 / math.sqrt(20)), (0.5 + 1 / math.sqrt(20)), 1]
    weights = [1.0 / 6, 5.0 / 6, 5.0 / 6, 1.0 / 6]
    n_points = len(q_points)

    u_point_val = numpy.zeros((3, u_n_inter * n_points))
    ud_point_val = numpy.zeros((3, u_n_inter * n_points))
    u_point_idx = numpy.zeros((u_n_inter * n_points, 1))
    q_u_point = numpy.zeros((u_n_inter * n_points, 1))

    n = 0
    for i in range(0, u_n_inter):
        us = u_knots[i + 2]
        uil = u_knots[i + 3] - u_knots[i + 2]
        my_list = range(n_points)
        for j in my_list:
            up = us + uil * q_points[j]
            q_u_point[n] = up  #
            u_base_vec, i_idx = spline_base_vec(u_knots, up, 0, True)
            u_base_vec_diff, i_idx = spline_base_vec(u_knots, up, 1, True)
            u_point_val[:, n] = u_base_vec
            ud_point_val[:, n] = u_base_vec_diff
            u_point_idx[n] = i_idx
            n += 1

    v_point_val = numpy.zeros((3, v_n_inter * n_points))
    vd_point_val = numpy.zeros((3, v_n_inter * n_points))
    v_point_idx = numpy.zeros((v_n_inter * n_points, 1))
    q_v_point = numpy.zeros((v_n_inter * n_points, 1))

    n = 0
    for i in range(v_n_inter):
        vs = v_knots[i + 2]
        vil = v_knots[i + 3] - v_knots[i + 2]
        my_list = range(n_points)
        for j in my_list:
            vp = vs + vil * q_points[j]
            q_u_point[n] = vp
            v_base_vec, j_idx = spline_base_vec(v_knots, vp, 0, True)
            v_base_vec_diff, j_idx = spline_base_vec(v_knots, vp, 1, True)
            v_point_val[:, n] = v_base_vec
            vd_point_val[:, n] = v_base_vec_diff
            v_point_idx[n] = j_idx
            n += 1

    # Matrix construction
    colv = numpy.zeros(9)
    data = numpy.zeros(9)
    data2 = numpy.zeros(9)
    row_m = numpy.zeros((v_n_inter * u_n_inter * n_points * n_points * 81))
    col_m = numpy.zeros((v_n_inter * u_n_inter * n_points * n_points * 81))
    data_m = numpy.zeros((v_n_inter * u_n_inter * n_points * n_points * 81))
    nnz_a = 0

    for i in range(v_n_inter):
        for k in range(n_points):
            v_point = v_point_val[:, i * n_points + k]
            vd_point = vd_point_val[:, i * n_points + k]
            j_idx = v_point_idx[i * n_points + k]
            for l in range(u_n_inter):
                for m in range(n_points):
                    u_point = u_point_val[:, l * n_points + m]
                    ud_point = ud_point_val[:, l * n_points + m]
                    i_idx = u_point_idx[l * n_points + m]
                    for n in range(0, 3):
                        # Hard-coded Kronecker product: vd = numpy.kron(vd_point, u_point)
                        data[3 * n:3 * (n + 1)] = vd_point[n] * u_point
                        # Hard-coded Kronecker product: ud = numpy.kron(v_point, ud_point)
                        data2[3 * n:3 * (n + 1)] = v_point[n] * ud_point
                        # column indices for data & data2
                        for p in range(0, 3):
                            colv[3 * n + p] = (j_idx + n) * u_n_basf + i_idx + p

                    # Hard-coded Outer product:
                    # Jacobian * weights[m] * weights[k] * (numpy.outer(ud, ud) + numpy.outer(vd, vd))
                    coef = weights[m] * weights[k] * compute_jacobian(q_u_point[l * n_points + m, 0],
                                                                      q_v_point[i * n_points + k, 0], a, b, c, d)
                    for n in range(0, 9):
                        row_m[nnz_a + 9 * n:nnz_a + 9 * (n + 1)] = colv
                        col_m[nnz_a + 9 * n:nnz_a + 9 * (n + 1)] = colv[n]
                        data_m[nnz_a + 9 * n:nnz_a + 9 * (n + 1)] = coef * (data[n] * data + data2[n] * data2)
                    nnz_a += 81

    mat_a = scipy.sparse.coo_matrix((data_m, (row_m, col_m)), shape=(u_n_basf * v_n_basf, u_n_basf * v_n_basf)).tocsr()

    return mat_a


def build_dense_reg_matrix(u_knots, v_knots):
    """
    Construction of the regularization matrix (A) to decrease variation of the terrain
    B z = b  ---> (B^T B + A)z = B^T b
    :param u_knots: vector of v-knots
    :param v_knots: vector of u-knots
    :return:
    """

    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3
    u_n_inter = len(u_knots) - 5
    v_n_inter = len(v_knots) - 5

    q_points = [0, (0.5 - 1 / math.sqrt(20)), (0.5 + 1 / math.sqrt(20)), 1]
    weights = [1.0 / 6, 5.0 / 6, 5.0 / 6, 1.0 / 6]
    n_points = len(q_points)

    u_point_val = numpy.zeros((u_n_basf, u_n_inter * n_points))
    ud_point_val = numpy.zeros((u_n_basf, u_n_inter * n_points))

    n = 0
    for i in range(0, u_n_inter):
        us = u_knots[i + 2]
        uil = u_knots[i + 3] - u_knots[i + 2]
        my_list = range(n_points)
        for j in my_list:
            up = us + uil * q_points[j]
            u_base_vec, i_idx = spline_base_vec(u_knots, up, 0, False)
            u_base_vec_diff, i_idx = spline_base_vec(u_knots, up, 1, False)
            u_point_val[:, n] = u_base_vec
            ud_point_val[:, n] = u_base_vec_diff
            n += 1

    v_point_val = numpy.zeros((v_n_basf, v_n_inter * n_points))
    vd_point_val = numpy.zeros((v_n_basf, v_n_inter * n_points))

    n = 0
    for i in range(v_n_inter):
        vs = v_knots[i + 2]
        vil = v_knots[i + 3] - v_knots[i + 2]
        my_list = range(n_points)
        for j in my_list:
            vp = vs + vil * q_points[j]
            v_base_vec, j_idx = spline_base_vec(v_knots, vp, 0, False)
            v_base_vec_diff, j_idx = spline_base_vec(v_knots, vp, 1, False)
            v_point_val[:, n] = v_base_vec
            vd_point_val[:, n] = v_base_vec_diff
            n += 1

    # Matrix construction
    mat_a = numpy.zeros((u_n_basf * v_n_basf, u_n_basf * v_n_basf))
    for i in range(v_n_inter):
        for k in range(n_points):
            v_point = v_point_val[:, i * n_points + k]
            vd_point = vd_point_val[:, i * n_points + k]
            for l in range(u_n_inter):
                for m in range(n_points):
                    u_point = u_point_val[:, l * n_points + m]
                    ud_point = ud_point_val[:, l * n_points + m]
                    vd = numpy.kron(vd_point, u_point)
                    ud = numpy.kron(v_point, ud_point)
                    mat_a += weights[m] * weights[k] * (numpy.outer(ud, ud) + numpy.outer(vd, vd))

    return mat_a


def build_reg_matrix(u_knots, v_knots, quad, sparse=True):
    """
    Construction of the regularization matrix (A) to decrease variation of the terrain
    B z = b  ---> (B^T B + A)z = B^T b
    :param u_knots: vector of v-knots
    :param v_knots: vector of u-knots
    :param quad: points defining quadrangle area (array)
    :param sparse: use sparse version or not
    :return:
    """

    if sparse is True:
        return build_sparse_reg_matrix(u_knots, v_knots, quad)
    else:
        return build_dense_reg_matrix(u_knots, v_knots)


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


def z_mat_to_bspline(u_knots, v_knots, z_mat, quad):
    """
    This function create B-Spline patch in raw format
    :param u_knots:
    :param v_knots:
    :param z_mat: iz = iu + nu * iv
    :param quad:
    """
    u_n_basf = len(u_knots) - 3
    v_n_basf = len(v_knots) - 3

    p0 = quad[:, 0]
    p1 = quad[:, 1]
    p2 = quad[:, 2]
    p3 = quad[:, 3]

    # Create list of poles from z_mat
    poles = [[[0.0, 0.0, 0.0] for i in range(v_n_basf)] for j in range(u_n_basf)]

    for j in range(0, v_n_basf):
        for i in range(0, u_n_basf):

            u = float(i) / (u_n_basf - 1)
            v = float(j) / (v_n_basf - 1)
            # Bi-linear interpolation of point?!
            xy = u * (v * p2 + (1 - v) * p1) + (1 - u) * (v * p3 + (1 - v) * p0)
            x_coord = xy[0, 0]
            y_coord = xy[1, 0]
            z_coord = z_mat[j * u_n_basf + i]
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

    return poles, tuple(u_knots), tuple(v_knots), tuple(u_mults), tuple(v_mults), u_deg, v_deg


def approx_svd(terrain_data, u_knots, v_knots, sparse=False, filter_thresh=0.001):
    """
    This function tries to approximate terrain data with B-Spline surface patches
    using SVD decomposition
    :param terrain_data: matrix of 3D terrain data
    :param u_knots: array of u knots
    :param v_knots: array of v knots
    :param sparse: use sparse matrices or not
    :param filter_thresh: threshold used in filtering S matrix
    :return: B-Spline patch
    """

    print('Assembling B matrix ...')
    start_time = time.time()
    mat_b, interval = build_ls_matrix(u_knots, v_knots, terrain_data, sparse)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    if sparse is True:
        mat_g = scipy.sparse.dok_matrix(terrain_data[:, 2])
    else:
        mat_g = numpy.matrix(terrain_data[:, 2])

    print('Computing SVD ...')
    start_time = time.time()
    if sparse is True:
        mat_u, mat_s, mat_v = numpy.linalg.svd(mat_b.todense(), full_matrices=False, compute_uv=True)
    else:
        mat_u, mat_s, mat_v = numpy.linalg.svd(mat_b, full_matrices=False, compute_uv=True)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Min(Matrix S): {0}, Max(Matrix S): {1}'.format(min(mat_s), max(mat_s)))

    print('Creating Si matrix ...')
    start_time = time.time()

    # Make computed matrices sparse to save memory
    if sparse is True:
        mat_u = scipy.sparse.dok_matrix(mat_u)
        mat_v = scipy.sparse.dok_matrix(mat_v)

    # Make compatible with Matlab code
    mat_v = mat_v.transpose()

    # Filter too small values and invert other values
    for key, value in enumerate(mat_s):
        if value < filter_thresh:
            mat_s[key] = 0.0
        else:
            mat_s[key] = 1 / value

    # rank = numpy.linalg.matrix_rank(mat_s)

    size = max(mat_s.shape)

    if sparse is True:
        mat_si = scipy.sparse.spdiags(mat_s, 0, size, size)
    else:
        mat_si = numpy.diagflat(mat_s)

    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Computing Z matrix ...')
    start_time = time.time()
    z_mat = mat_v * (mat_si.transpose() * (mat_u.transpose() * mat_g))
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    diff = eval_diff(mat_b, z_mat, mat_g)

    z_mat = z_mat_to_bspline(u_knots, v_knots, z_mat)
    return z_mat, diff


def approx_qr(terrain_data, u_knots, v_knots, sparse=False):
    """
    This function tries to approximate terrain data with B-Spline surface patches
    using QR decomposition
    :param terrain_data: matrix of 3D terrain data
    :param u_knots: array of u knots
    :param v_knots: array of v knots
    :param sparse:
    :return: B-Spline patch
    """

    # Own computation of approximation
    print('Creating B matrix ...')
    start_time = time.time()
    # b_mat, interval = build_ls_matrix(u_knots, v_knots, terrain_data, sparse)
    b_mat, interval = build_ls_matrix(u_knots, v_knots, terrain_data, True)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    g_mat = terrain_data[:, 2]

    # print('Computing QR ...')
    # start_time = time.time()
    # if sparse is True:
    #    q_mat, r_mat = numpy.linalg.qr(b_mat.todense(), mode='full')
    # else:
    #    q_mat, r_mat = numpy.linalg.qr(b_mat, mode='full')
    # end_time = time.time()
    # print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Computing Z matrix ...')
    start_time = time.time()
    # z_mat, diff = numpy.linalg.lstsq(b_mat.todense(), g_mat)[0]
    z_mat, diff = scipy.linalg.lstsq(b_mat.todense(), g_mat)
    # z_mat = numpy.linalg.lstsq(r_mat, q_mat.transpose() * g_mat)[0]
    # z_mat = scipy.linalg.solve_triangular(r_mat, q_mat.transpose() * g_mat)#[0]
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    # diff = eval_diff(b_mat, z_mat, g_mat)
    # diff = ((abs(g_mat - numpy.dot(b_mat, z_mat))).transpose()).tolist()[0]
    # z_mat_csr = (scipy.sparse.csr_matrix(z_mat))#.transpose()
    # diff = (b_mat.todense()* z_mat - g_mat).transpose().tolist()[0]
    # print(z_mat)
    # print(z_mat.shape)
    diff = diff.transpose().tolist()[0]
    z_mat = z_mat.transpose().tolist()[0]
    # diff = b_mat.dot(z_mat) - g_mat
    # diff = [val[0,0] for val in diff]

    z_mat = z_mat_to_bspline(u_knots, v_knots, z_mat)
    return z_mat, diff


def transform_points(quad, terrain_data):
    """
    Function computes corresponding (u,v) for (x,y)
    :param terrain_data: matrix of 3D terrain data
    :param quad: points defining quadrangle area (array)
    :return: XXX
    """
    mat_n = numpy.empty_like(quad)

    mat_n[:, 0] = quad[:, 0] - quad[:, 3]
    nt = mat_n[0, 0]
    mat_n[0, 0] = -mat_n[1, 0]
    mat_n[1, 0] = nt
    mat_n[:, 0] = mat_n[:, 0]/numpy.linalg.norm(mat_n[:, 0])

    for i in range(1, 4):
        mat_n[:, i] = quad[:, i] - quad[:, i-1]
        nt = mat_n[0, i]
        mat_n[0, i] = -mat_n[1, i]
        mat_n[1, i] = nt
        mat_n[:, i] = mat_n[:, i]/numpy.linalg.norm(mat_n[:, i])

    terrain_len = len(terrain_data)

    # Compute local coordinates and drop all points outside quadraangle
    param_terrain = numpy.empty_like(terrain_data)
    # indexing
    d0 = (terrain_data[:, 0:2] - quad[:, 0].transpose())
    d2 = (terrain_data[:, 0:2] - quad[:, 2].transpose())
    d3 = (terrain_data[:, 0:2] - quad[:, 3].transpose())
    d0_n0 = d0 * mat_n[:, 0]
    d0_n1 = d0 * mat_n[:, 1]
    d2_n2 = d2 * mat_n[:, 2]
    d3_n3 = d3 * mat_n[:, 3]
    u = numpy.divide(d0_n0, d0_n0 + d2_n2)
    v = numpy.divide(d0_n1, d0_n1 + d3_n3)

    h = -1
    for j in range(0,terrain_len):
        if (u[j] >= 0.0) and (u[j] <= 1.0) and (v[j] >= 0.0) and (v[j] <= 1.0):
            h += 1
            param_terrain[h, 0] = u[j]
            param_terrain[h, 1] = v[j]
            param_terrain[h, 2] = terrain_data[j, 2]

    param_terrain.resize(h+1, 3)

    uv = numpy.reshape(param_terrain[:, 0:2], 2*h+2).transpose()

    a = quad[:, 3]-quad[:, 2]
    b = quad[:, 0]-quad[:, 1]
    c = quad[:, 1]-quad[:, 2]
    d = quad[:, 0]-quad[:, 3]

    ldiag = numpy.zeros([2 * h + 1, 1])
    diag = numpy.zeros([2 * h + 2, 1])
    udiag = numpy.zeros([2 * h + 1, 1])

    # fixed point Newton iteration
    for i in range(0, 1):  # 1->5
        for j in range(0, h+1):
            mat_j = compute_jacobi(uv[2 * j, 0], uv[2 * j + 1, 0], a, b, c, d, -1)
            ldiag[2 * j, 0] = mat_j[1, 0]
            diag[2 * j, 0] = mat_j[0, 0]
            diag[2 * j + 1, 0] = mat_j[1, 1]
            udiag[2 * j, 0] = mat_j[0, 1]

        mat_jg = scipy.sparse.diags([ldiag[:, 0], diag[:, 0], udiag[:, 0]], [-1, 0, 1], format="csr")
        uv = uv - mat_jg.dot(uv)

    uv = uv.reshape([h + 1, 2])

    # Tresholding of the refined coordinates
    for j in range(0, h+1):
        if uv[j, 0] < 0:
            uv[:, 0] = 0
        elif uv[j, 0] > 1:
            uv[j, 0] = 1
        elif uv[j, 1] < 0:
            uv[:, 1] = 0
        elif uv[j, 1] > 1:
            uv[j, 1] = 1

    param_terrain[:,0:2] = uv

    return param_terrain


def compute_jacobi(u, v, a, b, c, d, output_type):
    """
    Compute Jacobian for local coordinates u & v
    :param u: value of parameter u
    :param v: value of parameter v
    :param a, b, c, d: vectors ( numpy array 2x1)
    :return:
    """

    mat_j = numpy.append(v * a + (1 - v) * b, u * c + (1 - u) * d, axis=1)

    if output_type is 1:
        return mat_j
    elif output_type is -1:
        return numpy.linalg.inv(mat_j)


def compute_jacobian(u, v, a, b, c, d):
    """
    Compute Jacobian for local coordinates u & v
    :param u: value of parameter u
    :param v: value of parameter v
    :param a, b, c, d: vectors ( numpy array 2x1)
    :return:
    """

    mat_j = numpy.append(v * a + (1 - v) * b, u * c + (1 - u) * d, axis=1)
    return numpy.linalg.det(mat_j)


def approx_chol(terrain_data, quad, u_knots, v_knots, sparse, filter_thresh):
    """
    This function tries to approximate terrain data with B-Spline surface patches
    using Cholesky decomposition
    :param terrain_data: matrix of 3D terrain data
    :param quad: points defining quadrangle area (array)
    :param u_knots: array of u knots
    :param v_knots: array of v knots
    :param sparse: if sparse matrix is used
    :param filter_thresh: threshold of filter
    :return: B-Spline patch
    """

    print('Transforming points to parametric space ...')
    start_time = time.time()
    param_terrain_data = transform_points(quad, terrain_data)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    # Own computation of approximation
    print('Creating B matrix ...')
    start_time = time.time()
    b_mat, interval = build_ls_matrix(u_knots, v_knots, param_terrain_data, sparse)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Creating A matrix ...')
    start_time = time.time()
    a_mat = build_reg_matrix(u_knots, v_knots, quad, sparse)
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    print('Computing B^T B matrix ...')
    start_time = time.time()
    if sparse is True:
        bb_mat = b_mat.transpose() * b_mat
    else:
        bb_mat = numpy.dot(b_mat.transpose(), b_mat)

    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    if sparse is True:
        # scipy.sparse.linalg.norm(bb_mat, 1)
        bb_norm = scipy.sparse.linalg.svds(bb_mat, k=1, ncv=10, tol=1e-4, which='LM', v0=None,
                                           maxiter=300, return_singular_vectors=False)
        # bb_mnsg = scipy.sparse.linalg.svds(bb_mat, k=1, ncv=10, tol=1e-2, which='SM', v0=None,
        # maxiter=300, return_singular_vectors=False)
        a_norm = scipy.sparse.linalg.svds(a_mat, k=1, ncv=10, tol=1e-4, which='LM', v0=None,
                                          maxiter=300, return_singular_vectors=False)
    else:
        bb_norm = numpy.linalg.norm(bb_mat)
        a_norm = numpy.linalg.norm(a_mat)

    r = filter_thresh
    c_mat = bb_mat + r * (bb_norm[0] / a_norm[0]) * a_mat

    g_vec = param_terrain_data[:, 2]
    b_vec = b_mat.transpose() * g_vec

    print('Computing Z coordinates ...')
    start_time = time.time()
    if sparse is True:
        z_vec = scipy.sparse.linalg.spsolve(c_mat, b_vec)
        diff = (numpy.matrix(b_mat * z_vec).transpose() - g_vec).tolist()[0]
    else:
        z_vec = scipy.linalg.solve(c_mat, b_vec)
        print(b_mat.shape)
        print(z_vec.shape)
        print((abs(g_vec - numpy.dot(b_mat, z_vec))).shape)
        diff = ((abs(g_vec - numpy.dot(b_mat, z_vec))).transpose()).tolist()[0]
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    z_vec = z_mat_to_bspline(u_knots, v_knots, z_vec, quad)

    return z_vec, diff


def eval_diff(b_mat, z_mat, g_mat):
    """
    # Compute difference between original terrain data and B-Spline surface
    :param b_mat: matrix of linear system
    :param z_mat: solution of the linear system
    :param g_mat: original RHS
    :return: list of differences
    """
    print('Computing differences ...')
    start_time = time.time()
    diff = ((abs(g_mat - numpy.dot(b_mat, z_mat))).transpose()).tolist()[0]
    end_time = time.time()
    print('Computed in {0} seconds.'.format(end_time - start_time))

    return diff


def approx(method, terrain_data, u_knots, v_knots, quad=None, sparse=False, conf=None):
    """
    This function tries to approximate terrain data with B-Spline surface patches
    :param method: method used for approximation
    :param terrain_data: matrix of 3D terrain data
    :param u_knots: array of u knots
    :param v_knots: array of v knots
    :param quad:
    :param sparse: sparse matrices will be used for computing
    :param conf: dictionary of other configuration specific for approximation method
    :return: B-Spline patch
    """

    conf = conf or {}

    if method == 'qr':
        return approx_qr(terrain_data, u_knots, v_knots, sparse)
    elif method == 'svd':
        return approx_svd(terrain_data, u_knots, v_knots, sparse, conf['threshold'])
    elif method == 'chol':
        return approx_chol(terrain_data, quad, u_knots, v_knots, sparse, conf['threshold'])
    else:
        raise TypeError("Wrong argument method: {0}".format(method))
