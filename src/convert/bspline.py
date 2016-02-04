"""
This module includes function for convertorting representation of
b-spline surface.
"""

import OCC
import OCC.gp
import OCC.TColgp
import OCC.TColStd
import OCC.Geom


class BSplineSurfConvertor(object):
    """
    Class used for converting SciPy representation
    of B-Spline surface to OpenCascade representation.
    """

    def __init__(self, tck):
        """
        Simple 'constructor' of BSplineSurfConvertor
        """
        super(BSplineSurfConvertor, self).__init__()
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
        weights = OCC.TColStd.TColStd_Array2OfReal(poles.LowerCol(), poles.ColLength(), poles.LowerRow(), poles.RowLength())
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
    convertor = BSplineSurfConvertor(scipy_bspline)
    occ_bspline = convertor.convert()
    return occ_bspline
