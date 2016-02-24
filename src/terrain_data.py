"""
This module is used for loading terrain data
"""

import yaml
import math
import sys
import colorsys

import OCC.BRep
import OCC.TopoDS
import OCC.BRepBuilderAPI
import OCC.BRepTools
import OCC.Prs3d
import OCC.Quantity
import OCC.Graphic3d
import OCC.Aspect
import OCC.gp

import convert.bspline


class TerrainData(object):
    """
    Class representing Terrain data (terrain, rivers, approximation, etc)
    """
    def __init__(self, yaml_file_name):
        """Constructor of TerrainData"""
        super(TerrainData, self).__init__()
        self.yaml_file_name = yaml_file_name
        self.conf = {}  # Configuration loaded from yaml file
        self.terrain_data = []  # Terrain data list of tuples (3d coordinates)
        self.point_count = 0
        self.tck = {}  # Patches of b-spline surface
        self.raw = {}
        self.tX = None  #
        self.tY = None  # Terrain X,Y,Z cache for numpy API
        self.tZ = None  #
        self.tW = None  # Visualization purpose
        self.min_x = -sys.maxsize
        self.max_x = sys.maxsize
        self.min_y = -sys.maxsize
        self.max_y = sys.maxsize
        self.size_x = 0
        self.size_y = 0
        self.diff_x = 0.0
        self.diff_y = 0.0
        self.dx = 0.0
        self.dy = 0.0
        self.grid = {}  # Cache for bilinear interpolation
        self.points = None
        self.rivers_data_2d = {}
        self.rivers_data_3d = {}
        self.area_borders_2d = {}
        self.area_borders_3d = {}
        self.bspline_surfaces = {}

    def load_conf_from_yaml(self):
        """
        Load configuration form yaml file
        """
        with open(self.yaml_file_name, 'r') as yaml_file:
            self.conf = yaml.load(yaml_file)

    def __post_process_terrain_data(self):
        """
        Try to post-process terrain data
        """
        self.point_count = len(self.terrain_data)
        self.tX = [item[0] for item in self.terrain_data]
        self.tY = [item[1] for item in self.terrain_data]
        self.tZ = [item[2] for item in self.terrain_data]

        def find_first_different(array):
            """
            Find index of first different item
            :param array: list of 3D points
            :return first index of different 3D point
            """
            a_index = 1
            first = array[0]
            while array[a_index] == first:
                a_index += 1
            a_index -= 1
            if a_index < len(array):
                return a_index
            else:
                return -1

        def find_first_same(array):
            """
            Find index of first same item
            :param array: list of 3D points
            :return index of first same 3D point
            """
            a_index = 1
            first = array[0]
            while array[a_index] != first:
                a_index += 1
            a_index -= 1
            if a_index < len(array):
                return a_index
            else:
                return -1

        self.min_x = min(self.tX)
        self.max_x = max(self.tX)
        self.min_y = min(self.tY)
        self.max_y = max(self.tY)
        # Try to compute size of 2d array according repeated values
        if self.tX[0] == self.tX[1]:
            self.size_x = find_first_different(self.tX) + 1
        else:
            self.size_x = find_first_same(self.tX) + 1
        self.size_y = self.point_count / self.size_x
        self.diff_x = self.max_x - self.min_x
        self.diff_y = self.max_y - self.min_y
        self.dx = self.diff_x / float(self.size_x - 1)
        self.dy = self.diff_y / float(self.size_y - 1)
        self.grid = {}
        # Add terrain points to grid
        self.points = [(0.0, 0.0, 0.0) for i in range(0, self.point_count)]
        for index in range(0, self.point_count):
            i = int(math.floor((self.terrain_data[index][0] - self.min_x) / self.dx))
            j = int(math.floor((self.terrain_data[index][1] - self.min_y) / self.dy))
            x_coord = self.terrain_data[index][0]
            y_coord = self.terrain_data[index][1]
            z_coord = self.terrain_data[index][2]
            self.grid[(i, j)] = (x_coord, y_coord, z_coord)
            # Transform x, y coordinates to range <0, 1>
            x_coord = (x_coord - self.min_x) / self.diff_x
            y_coord = (y_coord - self.min_y) / self.diff_y
            self.points[index] = (x_coord, y_coord, z_coord)

    def load_terrain(self):
        """
        Try to load data of terrain
        """
        with open(self.conf['terrain'], 'r') as data_file:
            for line in data_file:
                self.terrain_data.append(tuple(float(item) for item in line.split()))
        self.__post_process_terrain_data()

    def load_rivers(self):
        """
        Try to load data of rivers
        """
        with open(self.conf['rivers'], 'r') as data_file:
            for line in data_file:
                items = line.split()
                river_id = int(items[-1])
                river_coord = (float(items[0]), float(items[1]))
                try:
                    river = self.rivers_data_2d[river_id]
                except KeyError:
                    river = self.rivers_data_2d[river_id] = []
                river.append(river_coord)

    def load_area(self):
        """
        Try to load definition of area border
        """
        with open(self.conf['area'], 'r') as area_file:
            self.area_borders_2d[0] = []
            for line in area_file:
                items = line.split()
                # Only one border ATM
                self.area_borders_2d[0].append(tuple(float(item) for item in items[1:]))

    def approximate_terrain(self):
        """
        Try to approximate terrain with bspline surface
        """
        if self.conf['approximation']['solver'] == 'scipy':
            from scipy import interpolate
            tck, fp, ior, msg = interpolate.bisplrep(self.tX, self.tY, self.tZ, kx=5, ky=5, full_output=1)
            self.tck[(self.min_x, self.min_y, self.max_x, self.max_y)] = tck
            # Compute difference between original terrain data and B-Spline surface
            self.tW = [abs(it[2] - interpolate.bisplev(it[0], it[1], tck)) for it in self.terrain_data]
        elif self.conf['approximation']['solver'] == 'raw':
            import approx.terrain
            import numpy
            u_knots = approx.terrain.gen_knots(self.conf['approximation']['u_knots_num'])
            v_knots = approx.terrain.gen_knots(self.conf['approximation']['v_knots_num'])
            terrain = numpy.matrix(self.points)
            poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = approx.terrain.approx(terrain, u_knots, v_knots)
            # Transform x, y coordinates of poles back to original range,
            # because x, y coordinates were transformed to range <0, 1>
            for i in range(0, len(poles)):
                for j in range(0, len(poles[0])):
                    x_coord = self.min_x + self.diff_x * poles[i][j][0]
                    y_coord = self.min_y + self.diff_y * poles[i][j][1]
                    poles[i][j] = (x_coord, y_coord, poles[i][j][2])
            raw = (poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
            self.raw[(self.min_x, self.min_y, self.max_x, self.max_y)] = raw
            # TODO: Compute difference between original terrain data and B-Spline surface
            self.tW = [0.0 for it in self.terrain_data]

    def approximate_2d_borders(self):
        """
        Try to approximate z coordinates of borders using terrain data
        """

        for border_id, border_2d in self.area_borders_2d.items():
            self.area_borders_3d[border_id] = []
            for border_pnt in border_2d:
                # Compute indexes to the grid first
                i_index = int(math.floor((border_pnt[0] - self.min_x) / self.dx))
                j_index = int(math.floor((border_pnt[1] - self.min_y) / self.dy))
                # Compute weights for bilineral interpolation
                kx_weight = (border_pnt[0] - (self.min_x + self.dx * i_index)) / self.dx
                ky_weight = (border_pnt[1] - (self.min_y + self.dy * j_index)) / self.dy
                z1_coord = self.grid[(i_index, j_index)][2]
                z2_coord = self.grid[(i_index + 1, j_index)][2]
                z3_coord = self.grid[(i_index, j_index + 1)][2]
                z4_coord = self.grid[(i_index + 1, j_index + 1)][2]
                z12_coord = (1.0 - kx_weight) * z1_coord + kx_weight * z2_coord
                z34_coord = (1.0 - kx_weight) * z3_coord + kx_weight * z4_coord
                z_coord = (1.0 - ky_weight) * z12_coord + ky_weight * z34_coord
                self.area_borders_3d[border_id].append((border_pnt[0], border_pnt[1], z_coord))

    def approximate_2d_rivers(self):
        """
        Try to approximate z coordinates of rivers using terrain data
        """

        for river_id, river2d in self.rivers_data_2d.items():
            self.rivers_data_3d[river_id] = []
            last_z = sys.maxsize
            for river_pnt in river2d:
                # Compute indexes to the grid first
                i_index = int(math.floor((river_pnt[0] - self.min_x) / self.dx))
                j_index = int(math.floor((river_pnt[1] - self.min_y) / self.dy))
                # Compute weights for bilineral interpolation
                kx_weight = (river_pnt[0] - (self.min_x + self.dx * i_index)) / self.dx
                ky_weight = (river_pnt[1] - (self.min_y + self.dy * j_index)) / self.dy
                z1_coord = self.grid[(i_index, j_index)][2]
                z2_coord = self.grid[(i_index + 1, j_index)][2]
                z3_coord = self.grid[(i_index, j_index + 1)][2]
                z4_coord = self.grid[(i_index + 1, j_index + 1)][2]
                z12_coord = (1.0 - kx_weight) * z1_coord + kx_weight * z2_coord
                z34_coord = (1.0 - kx_weight) * z3_coord + kx_weight * z4_coord
                z_coord = (1.0 - ky_weight) * z12_coord + ky_weight * z34_coord
                # Cut out too big z values
                if last_z < z_coord:
                    z_coord = last_z
                last_z = z_coord
                self.rivers_data_3d[river_id].append((river_pnt[0], river_pnt[1], z_coord))

    @staticmethod
    def create_ogl_group(display):
        """
        Create a group that will store an OpenGL buffer
        :param display: OCC display
        :return tuple of OCC 3D presentation and current OCC 3D group
        """
        presentation = OCC.Prs3d.Prs3d_Presentation(display._struc_mgr)
        group = OCC.Prs3d.Prs3d_Root_CurrentGroup(presentation.GetHandle()).GetObject()
        return presentation, group

    def display_terrain(self):
        """
        Try to display terrain
        """

        # Initialize displaying
        from OCC.Display.SimpleGui import init_display
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.EraseAll()
        # Draw all bspline surfaces
        for bspline in self.bspline_surfaces.values():
            display.DisplayShape(bspline.GetHandle(), update=True)

        # Draw points of terrain
        a_presentation, group = self.create_ogl_group(display)
        black = OCC.Quantity.Quantity_Color(OCC.Quantity.Quantity_NOC_BLACK)
        asp = OCC.Graphic3d.Graphic3d_AspectLine3d(black, OCC.Aspect.Aspect_TOL_SOLID, 1)

        gg = OCC.Graphic3d.Graphic3d_ArrayOfPoints(self.point_count,
                                                   True,  # hasVColors
                                                   )

        max_diff = max(self.tW)
        idx = 1
        for point in self.terrain_data:
            pnt = OCC.gp.gp_Pnt(point[0], point[1], point[2])
            gg.AddVertex(pnt)
            # create the point, with a random color
            if max_diff > 0.0:
                diff = self.tW[idx - 1] / max_diff
            else:
                diff = 0.0
            rgb = colorsys.hsv_to_rgb(diff, 1.0, 1.0)
            gg.SetVertexColor(idx, rgb[0], rgb[1], rgb[2])
            idx += 1

        group.SetPrimitivesAspect(asp.GetHandle())
        group.AddPrimitiveArray(gg.GetHandle())
        a_presentation.Display()

        start_display()

    def output_approx_data(self):
        """
        Try to output approximated data to BREP file format
        """

        sewing = OCC.BRepBuilderAPI.BRepBuilderAPI_Sewing(0.01, True, True, True, False)
        sewing.SetFloatingEdgesMode(True)

        error = 1e-6

        for key, scipy_bspline in self.tck.items():
            occ_bspline = convert.bspline.scipy_to_occ(scipy_bspline)
            self.bspline_surfaces[key] = occ_bspline
            face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(occ_bspline.GetHandle(), error).Shape()
            sewing.Add(face)

        for key, raw_bspline in self.raw.items():
            poles = raw_bspline[0]
            u_knots = raw_bspline[1]
            v_knots = raw_bspline[2]
            u_mults = raw_bspline[3]
            v_mults = raw_bspline[4]
            u_deg = raw_bspline[5]
            v_deg = raw_bspline[6]
            occ_bspline = convert.bspline.raw_to_occ(poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
            self.bspline_surfaces[key] = occ_bspline
            face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(occ_bspline.GetHandle(), error).Shape()
            sewing.Add(face)

        sewing.Perform()
        sewing_shape = sewing.SewedShape()

        # shell = OCC.TopoDS.topods_Shell(sewing_shape)

        # make_solid = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeSolid()
        # make_solid.Add(shell)

        # solid = make_solid.Solid()

        # builder.MakeSolid(solid)
        # builder.Add(solid, shell)

        # compound = TopoDS_Compound()

        # builder = OCC.BRep.BRep_Builder()
        # builder.MakeCompound(compound)
        # builder.Add(compound, solid)

        # OCC.BRepTools.breptools_Write(compound, self.conf['output'])

        OCC.BRepTools.breptools_Write(sewing_shape, self.conf['output'])
