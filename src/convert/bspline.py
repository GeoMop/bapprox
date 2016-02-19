"""
This module includes function for converting representation of
b-spline surface.
"""

import OCC
import OCC.gp
import OCC.TColgp
import OCC.TColStd
import OCC.Geom


class BSplineSurfConvertorSciPy2OCC(object):
    """
    Class used for converting SciPy representation
    of B-Spline surface to OpenCascade representation.
    """

    def __init__(self, tck):
        """
        Simple 'constructor' of BSplineSurfConvertor
        """
        super(BSplineSurfConvertorSciPy2OCC, self).__init__()
        self.tck = tck

    @staticmethod
    def __unique_values(t_values):
        """
        This computes number of unique T values
        """
        count = 0
        last_val = None
        for t_val in t_values:
            if last_val is None or last_val != t_val:
                count += 1
            last_val = t_val
        return count

    @staticmethod
    def __generate_poles(tx_values, ty_values, cont, col_len, row_len):
        """
        This function generates OCC 2D array of poles
        """
        u_indexes = range(0, col_len)
        v_indexes = range(0, row_len)

        min_x = min(tx_values)
        max_x = max(tx_values)
        min_y = min(ty_values)
        max_y = max(ty_values)
        diff_x = float(max_x - min_x) / (col_len - 1)
        diff_y = float(max_y - min_y) / (row_len - 1)

        poles = OCC.TColgp.TColgp_Array2OfPnt(1, col_len, 1, row_len)
        # Set poles of b-spline surface
        c_i = 0
        for key_u in u_indexes:
            for key_v in v_indexes:
                x_coord = min_x + diff_x * key_u
                y_coord = min_y + diff_y * key_v
                z_coord = cont[c_i]
                point = OCC.gp.gp_Pnt(x_coord, y_coord, z_coord)
                poles.SetValue(key_u + 1, key_v + 1, point)
                c_i += 1
        return poles

    @staticmethod
    def __generate_weights(poles, weight_fnc=lambda u, v: 1.0):
        """
        This function generates array of weights for poles, based on u and v.
        When no weight function is defined, then it generates only default values (1.0)
        """
        weights = OCC.TColStd.TColStd_Array2OfReal(poles.LowerCol(), poles.ColLength(),
                                                   poles.LowerRow(), poles.RowLength())
        for key_u in range(poles.LowerCol(), poles.ColLength()+1):
            for key_v in range(poles.LowerRow(), poles.RowLength()+1):
                weights.SetValue(key_u, key_v, weight_fnc(key_u, key_v))
        return weights

    @staticmethod
    def __generate_knots(t_values, knots_len):
        """
        This function generates OCC 1D array of knots from T-values
        """
        knots = OCC.TColStd.TColStd_Array1OfReal(1, knots_len)
        last_val = None
        key_i = 1
        for t_val in t_values:
            value = (t_val - t_values[0]) / (t_values[-1] - t_values[0])
            if last_val is None or last_val != value:
                knots.SetValue(key_i, value)
                key_i += 1
            last_val = value
        return knots

    @staticmethod
    def __generate_mults(t_values, mult_len):
        """
        This function generates OCC 1D array of multiplicities from T-values
        """
        mults = OCC.TColStd.TColStd_Array1OfInteger(1, mult_len)
        mult_i = 1
        last_val = None
        mult = 0
        for t_val in t_values:
            if last_val is None:
                mult = 1
            elif last_val == t_val:
                mult += 1
            else:
                mults.SetValue(mult_i, mult)
                mult_i += 1
                mult = 1
            last_val = t_val
        mults.SetValue(mult_i, mult)
        return mults

    def convert(self):
        """
        This method convert do the conversion.
        """
        tx_values = self.tck[0]
        ty_values = self.tck[1]
        cont = self.tck[2]
        udeg = self.tck[3]
        vdeg = self.tck[4]

        # Try to convert SciPY B-Spline description to OCC B-Spline description
        col_len = len(tx_values) - udeg - 1
        row_len = len(ty_values) - vdeg - 1

        poles = self.__generate_poles(tx_values, ty_values, cont, col_len, row_len)

        weights = self.__generate_weights(poles)

        # This have to hold
        uknot_len = umult_len = self.__unique_values(tx_values)
        vknot_len = vmult_len = self.__unique_values(ty_values)

        # Set knots of b-spline surface
        uknots = self.__generate_knots(tx_values, uknot_len)
        vknots = self.__generate_knots(ty_values, vknot_len)

        # Set multiplicities of b-spline surface
        umult = self.__generate_mults(tx_values, umult_len)
        vmult = self.__generate_mults(ty_values, vmult_len)

        return OCC.Geom.Geom_BSplineSurface(poles, weights, uknots, vknots, umult, vmult, udeg, vdeg, False, False)


def scipy_to_occ(scipy_bspline):
    """
    This function converts SciPy representation of B-Spline surface to
    OCC representation.
    """
    convertor = BSplineSurfConvertorSciPy2OCC(scipy_bspline)
    occ_bspline = convertor.convert()
    return occ_bspline


class BSplineSurfConvertorRaw2OCC(object):
    """
    Class used for converting raw representation
    of B-Spline surface (tuple of poles, knots, multiplicities
    and degrees) to OpenCascade representation.
    """

    def __init__(self, poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg):
        """
        Simple 'constructor' of BSplineSurfConvertorRaw2OCC
        """
        super(BSplineSurfConvertorRaw2OCC, self).__init__()
        self.poles = poles
        self.u_knots = u_knots
        self.v_knots = v_knots
        self.u_mults = u_mults
        self.v_mults = v_mults
        self.u_deg = u_deg
        self.v_deg = v_deg

    def __generate_poles(self):
        """
        This method generates OCC array of points
        """
        col_len = len(self.poles)
        row_len = len(self.poles[0])
        poles = OCC.TColgp.TColgp_Array2OfPnt(1, col_len, 1, row_len)
        for key_u in range(0, col_len):
            for key_v in range(0, row_len):
                # print(self.poles[key_u][key_v])
                x_coord = self.poles[key_u][key_v][0]
                y_coord = self.poles[key_u][key_v][1]
                z_coord = self.poles[key_u][key_v][2]
                # print(x_coord, y_coord, z_coord)
                point = OCC.gp.gp_Pnt(x_coord, y_coord, z_coord)
                poles.SetValue(key_u + 1, key_v + 1, point)
        return poles

    def __generate_knots(self):
        """
        This function return tuple of two OCC knot arrays (u_knots, v_knots)
        :return: tuple of two OCC knot arrays
        """
        u_knots_len = len(self.u_knots)
        v_knots_len = len(self.v_knots)
        u_knots = OCC.TColStd.TColStd_Array1OfReal(1, u_knots_len)
        v_knots = OCC.TColStd.TColStd_Array1OfReal(1, v_knots_len)
        for key_i in range(0, u_knots_len):
            u_knots.SetValue(key_i + 1, self.u_knots[key_i])
        for key_i in range(0, v_knots_len):
            v_knots.SetValue(key_i + 1, self.v_knots[key_i])
        return u_knots, v_knots

    def __generate_mults(self):
        """
        Thus function returns tuple of two OCC multiplicities array (u_mults, v_mults)
        :return:
        """
        u_mults_len = len(self.u_mults)
        v_mults_len = len(self.v_mults)
        u_mults = OCC.TColStd.TColStd_Array1OfInteger(1, u_mults_len)
        v_mults = OCC.TColStd.TColStd_Array1OfInteger(1, v_mults_len)
        for key_i in range(0, u_mults_len):
            u_mults.SetValue(key_i + 1, self.u_mults[key_i])
        for key_i in range(0, v_mults_len):
            v_mults.SetValue(key_i + 1, self.v_mults[key_i])
        return u_mults, v_mults

    def convert(self):
        """
        This method do the conversion
        """
        poles = self.__generate_poles()
        u_knots, v_knots = self.__generate_knots()
        u_mults, v_mults = self.__generate_mults()
        return OCC.Geom.Geom_BSplineSurface(poles, u_knots, v_knots, u_mults, v_mults, self.u_deg, self.v_deg, False, False)


def raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg):
    """
    This function converts raw representation of B-Spline surface to
    OCC representation.
    :param poles: Iterable param containing poles
    :param u_knots: Iterable param containing
    :param v_knots:
    :param u_mults:
    :param v_mults:
    :param u_deg:
    :param v_deg:
    :return:
    """
    convertor = BSplineSurfConvertorRaw2OCC(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
    occ_bspline = convertor.convert()
    return occ_bspline
