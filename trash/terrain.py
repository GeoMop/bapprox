"""
Trash code
"""

def build_reg_matrix(u_knots, v_knots, quad, sparse=True):
    """

    :param u_knots:
    :param v_knots:
    :param quad:
    :param sparse:
    :return:
    """
    # u_point_valx = numpy.zeros((u_n_basf, u_n_inter * n_points))
    # ud_point_valx = numpy.zeros((u_n_basf, u_n_inter * n_points))
    #
    # n = 0
    # for i in range(0, u_n_inter):
    #     us = u_knots[i+2]
    #     uil = u_knots[i+3] - u_knots[i+2]
    #     my_list = range(n_points)
    #     for j in my_list:
    #         up = us + uil * q_points[j]
    #         u_base_vec, i_idx  = spline_base_vec(u_knots, up, 0, 'false')
    #         u_base_vec_diff, i_idx  = spline_base_vec(u_knots, up, 1, 'false')
    #         u_point_valx[:, n] = u_base_vec
    #         ud_point_valx[:, n] = u_base_vec_diff
    #         n += 1
    #
    # v_point_valx = numpy.zeros((v_n_basf, v_n_inter * n_points))
    # vd_point_valx = numpy.zeros((v_n_basf, v_n_inter * n_points))
    #
    # n = 0
    # for i in range(v_n_inter):
    #     vs = v_knots[i+2]
    #     vil = v_knots[i+3] - v_knots[i+2]
    #     my_list = range(n_points)
    #     for j in my_list:
    #         vp = vs + vil * q_points[j]
    #         v_base_vec, j_idx  = spline_base_vec(v_knots, vp, 0, 'false')
    #         v_base_vec_diff, j_idx  = spline_base_vec(v_knots, vp, 1, 'false')
    #         v_point_valx[:, n] = v_base_vec
    #         vd_point_valx[:, n] = v_base_vec_diff
    #         n += 1
    #
    #     # Matrix construction
    # mat_x = numpy.zeros((u_n_basf * v_n_basf, u_n_basf * v_n_basf))
    # for i in range(v_n_inter):
    #     for k in range(n_points):
    #         v_point = v_point_valx[:, i * n_points + k]
    #         vd_point = vd_point_valx[:, i * n_points + k]
    #         for l in range(u_n_inter):
    #             for m in range(n_points):
    #                 u_point = u_point_valx[:, l * n_points + m]
    #                 ud_point = ud_point_valx[:, l * n_points + m]
    #                 vd = numpy.kron(vd_point, u_point)
    #                 ud = numpy.kron(v_point, ud_point)
    #                 mat_x += weights[m] * weights[k] * (numpy.outer(ud, ud) + numpy.outer(vd, vd))
    #
    # print(numpy.linalg.norm(mat_a.todense() - mat_x))