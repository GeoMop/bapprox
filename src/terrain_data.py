"""
This module is used for loading terrain data
"""

import yaml
import math
import sys
import colorsys
import time

import OCC.BRep
import OCC.TopoDS
import OCC.BRepBuilderAPI
import OCC.BRepTools
import OCC.BRepAlgoAPI
import OCC.Prs3d
import OCC.Quantity
import OCC.Graphic3d
import OCC.Aspect
import OCC.gp
import OCC.BRepLib
import OCC.GeomAPI

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
        self.max_diff = 0.0
        self.min_x = -sys.maxsize
        self.max_x = sys.maxsize
        self.min_y = -sys.maxsize
        self.max_y = sys.maxsize
        self.min_z = -sys.maxsize
        self.max_z = sys.maxsize
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
        self.rivers_curves_3d = {}
        self.area_borders_2d = {}
        self.area_borders_3d = {}
        self.area_volumes = {}
        self.bspline_surfaces = {}
        self.shell = None
        self.volume = None
        self.sewing_shape = None
        self.fractures_shape = None
        self.quad_x = []
        self.quad_y = []


    def load_conf_from_yaml(self):
        """
        Load configuration form yaml file
        """
        # Yaml file is loaded to the dictionary in two lines of code :-)
        with open(self.yaml_file_name, 'r') as yaml_file:
            self.conf = yaml.load(yaml_file)

        # When some values are not set, then set default values
        # Terrain
        try:
            self.conf['terrain']
        except KeyError:
            self.conf['terrain'] = {}
        try:
            self.conf['terrain']['extrude_diff']
        except KeyError:
            self.conf['terrain']['extrude_diff'] = 0.0
        try:
            self.conf['terrain']['approximation']
        except KeyError:
            self.conf['terrain']['approximation'] = {}
        try:
            self.conf['terrain']['approximation']['solver']
        except KeyError:
            self.conf['terrain']['approximation']['solver'] = {}
        try:
            self.conf['terrain']['approximation']['solver']['method']
        except KeyError:
            self.conf['terrain']['approximation']['solver']['method'] = 'scipy'
        try:
            self.conf['terrain']['approximation']['solver']['sparse']
        except KeyError:
            self.conf['terrain']['approximation']['solver']['sparse'] = True
        try:
            self.conf['terrain']['approximation']['u_knots_num']
        except KeyError:
            self.conf['terrain']['approximation']['u_knots_num'] = 15
        try:
            self.conf['terrain']['approximation']['v_knots_num']
        except KeyError:
            self.conf['terrain']['approximation']['v_knots_num'] = 15
        try:
            self.conf['terrain']['approximation']['differences']
        except KeyError:
            self.conf['terrain']['approximation']['differences'] = False
        try:
            self.conf['terrain']['approximation']['output_differences']
        except KeyError:
            self.conf['terrain']['approximation']['output_differences'] = None
        try:
            self.conf['terrain']['approximation']['quad_x']
        except KeyError:
            self.conf['terrain']['approximation']['quad_x'] = None
        try:
            self.conf['terrain']['approximation']['quad_y']
        except KeyError:
            self.conf['terrain']['approximation']['quad_y'] = None
        # Output
        try:
            self.conf['output']
        except KeyError:
            self.conf['output'] = 'terrain.brep'
        # Area
        try:
            self.conf['area']
        except KeyError:
            self.conf['area'] = {}
        # Rivers
        try:
            self.conf['rivers']
        except KeyError:
            self.conf['rivers'] = {}
        # Fractures
        try:
            self.conf['fractures']
        except KeyError:
            self.conf['fractures'] = {}
        # Display results
        try:
            self.conf['display']
        except KeyError:
            self.conf['display'] = {}
        try:
            self.conf['display']['surface']
        except KeyError:
            self.conf['display']['surface'] = True
        try:
            self.conf['display']['terrain']
        except KeyError:
            self.conf['display']['terrain'] = False
        try:
            self.conf['display']['rivers']
        except KeyError:
            self.conf['display']['rivers'] = True
        try:
            self.conf['display']['fractures']
        except KeyError:
            self.conf['display']['fractures'] = False

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
        self.min_z = min(self.tZ)
        self.max_z = max(self.tZ)
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
            #x_coord = (x_coord - self.min_x) / self.diff_x
            #y_coord = (y_coord - self.min_y) / self.diff_y
            self.points[index] = (x_coord, y_coord, z_coord)
            #print(self.min_x,self.max_x,self.min_y,self.max_y)

    def load_terrain(self):
        """
        Try to load data of terrain
        """
        with open(self.conf['terrain']['input'], 'r') as data_file:
            for line in data_file:
                self.terrain_data.append(tuple(float(item) for item in line.split()))
        self.__post_process_terrain_data()

    def load_rivers(self):
        """
        Try to load data of rivers
        """
        if 'rivers' in self.conf and 'input' in self.conf['rivers']:
            with open(self.conf['rivers']['input'], 'r') as data_file:
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
        if 'area' in self.conf and 'input' in self.conf['area']:
            with open(self.conf['area']['input'], 'r') as area_file:
                self.area_borders_2d[0] = []
                for line in area_file:
                    items = line.split()
                    # Only one border ATM
                    self.area_borders_2d[0].append(tuple(float(item) for item in items[1:]))

    def load_fractures(self):
        """
        Try to load BREP file with fractures.
        :return:
        """
        if 'fractures' in self.conf and 'input' in self.conf['fractures']:
            self.fractures_shape = OCC.TopoDS.TopoDS_Shape()
            builder = OCC.BRep.BRep_Builder()
            OCC.BRepTools.breptools_Read(self.fractures_shape, self.conf['fractures']['input'], builder)

    def approximate_terrain(self):
        """
        Try to approximate terrain with bspline surface
        """
        solver_method = self.conf['terrain']['approximation']['solver']['method']
        sparse = self.conf['terrain']['approximation']['solver']['sparse']
        threshold = self.conf['terrain']['approximation']['solver']['threshold']
        u_knots_num = self.conf['terrain']['approximation']['u_knots_num']
        v_knots_num = self.conf['terrain']['approximation']['v_knots_num']
        comp_diffs = self.conf['terrain']['approximation']['differences']
        output_diff = self.conf['terrain']['approximation']['output_differences']
        quad_x = self.conf['terrain']['approximation']['quad_x']
        quad_y = self.conf['terrain']['approximation']['quad_y']

        if solver_method == 'scipy':
            from scipy import interpolate
            print('SciPy approximation ...')
            start_time = time.time()
            if u_knots_num == 0 and v_knots_num == 0:
                tck, fp, ior, msg = interpolate.bisplrep(self.tX, self.tY, self.tZ, kx=5, ky=5, full_output=1)
            else:
                tck, fp, ior, msg = interpolate.bisplrep(self.tX, self.tY, self.tZ, kx=4, ky=4,
                                                         nxest=u_knots_num,
                                                         nyest=v_knots_num,
                                                         full_output=1)
            end_time = time.time()
            print('Computed in {0} seconds with WSoSR: {1}.'.format(end_time - start_time, fp))
            if ior > 0:
                print('Warning({0}): {1}.'.format(ior, msg))
            self.tck[(self.min_x, self.min_y, self.max_x, self.max_y)] = tck
            # Compute difference between original terrain data and B-Spline surface
            if comp_diffs is True:
                print('Computing differences ...')
                start_time = time.time()
                self.tW = [abs(it[2] - interpolate.bisplev(it[0], it[1], tck)) for it in self.terrain_data]
                end_time = time.time()
                print('Computed in {0} seconds.'.format(end_time - start_time))
        elif solver_method in ['qr', 'svd', 'chol']:
            import approx.terrain
            import numpy
            u_knots = approx.terrain.gen_knots(u_knots_num)
            v_knots = approx.terrain.gen_knots(v_knots_num)
            terrain = numpy.matrix(self.points)
            # FIXME: When quad area is not defined, then define quad from max na min values of x coordinates
            # NOTE: This is not effective
            if quad_x is None or quad_y is None:
                quad = numpy.matrix([[self.min_x, self.max_x, self.max_x, self.min_x],
                                     [self.min_y, self.min_y, self.max_y, self.max_y]])
            else:
                quad = numpy.matrix([quad_x, quad_y])

            # Do own B-Spline approximation o terrain data
            if solver_method == 'chol':
                raw, diffs = approx.terrain.approx(solver_method, terrain, u_knots, v_knots, quad, sparse,
                                                   {'threshold': threshold})
            elif solver_method == 'svd':
                raw, diffs = approx.terrain.approx(solver_method, terrain, u_knots, v_knots, sparse=sparse,
                                                   conf={'threshold': threshold})
            else:
                raw, diffs = approx.terrain.approx(solver_method, terrain, u_knots, v_knots, sparse=sparse)
            poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg = raw
            if comp_diffs is True:
                self.tW = diffs
            # Transform x, y coordinates of poles back to original range,
            # because x, y coordinates were transformed to range <0, 1>
            if solver_method != 'chol':
                for i in range(0, len(poles)):
                   for j in range(0, len(poles[0])):
                       x_coord = self.min_x + self.diff_x * poles[i][j][0]
                       y_coord = self.min_y + self.diff_y * poles[i][j][1]
                       poles[i][j] = (x_coord, y_coord, poles[i][j][2])
                raw = (poles, u_knots, v_knots, u_mults, v_mults, u_deg, v_deg)
            self.raw[(self.min_x, self.min_y, self.max_x, self.max_y)] = raw

        if comp_diffs is True:
            #self.max_diff = self.tW.max()
            self.max_diff = max(self.tW)
            print('Max difference {0}'.format(self.max_diff))
            if output_diff is not None:
                self.output_diffs_to_csv_file(output_diff)

    def output_diffs_to_csv_file(self, output_diff_file):
        """
        :return None
        """
        with open(output_diff_file, 'w') as csv_diff_file:
            for idx, diff_val in enumerate(self.tW):
                x_coord, y_coord, z_coord = self.terrain_data[idx]
                csv_diff_file.write("{0} {1} {2} {3}\n".format(str(x_coord), str(y_coord), str(z_coord), str(diff_val)))

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
                # Compute weights for bi-linear interpolation
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

    def create_volume_from_area(self, extrude_diff):
        """
        Try to create volume from shape defined in area data file
        :param extrude_diff: distance between upper and lower cap of volume
        :return:
        """

        for border_id, border_2d in self.area_borders_2d.items():

            sewing = OCC.BRepBuilderAPI.BRepBuilderAPI_Sewing(0.01, True, True, True, False)
            sewing.SetFloatingEdgesMode(True)

            min_z = self.min_z + extrude_diff - 10.0
            max_z = self.max_z + 10.0

            cap1_wire = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeWire()
            cap2_wire = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeWire()
            prev_point = None
            for pnt_id, border_pnt in enumerate(border_2d):
                if prev_point is not None:
                    point0 = prev_point
                    point = point1 = OCC.gp.gp_Pnt(border_pnt[0], border_pnt[1], min_z)
                    point2 = OCC.gp.gp_Pnt(border_pnt[0], border_pnt[1], max_z)
                    point3 = OCC.gp.gp_Pnt(prev_point.X(), prev_point.Y(), max_z)
                    # Create wire (border of one face)
                    wire = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeWire()
                    # 0
                    edge0 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point0, point1)
                    wire.Add(edge0.Edge())
                    cap1_wire.Add(edge0.Edge())
                    # 1
                    edge1 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point1, point2)
                    wire.Add(edge1.Edge())
                    # 2
                    edge2 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point2, point3)
                    wire.Add(edge2.Edge())
                    cap2_wire.Add(edge2.Edge())
                    # 3
                    edge3 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point3, point0)
                    wire.Add(edge3.Edge())
                    # Create face from the wire
                    face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
                    sewing.Add(face.Shape())
                else:
                    first_point = point = OCC.gp.gp_Pnt(border_pnt[0], border_pnt[1], min_z)
                prev_point = point

            # Connect last point with first point
            point0 = OCC.gp.gp_Pnt(border_pnt[0], border_pnt[1], min_z)
            point1 = first_point
            point2 = OCC.gp.gp_Pnt(first_point.X(), first_point.Y(), max_z)
            point3 = OCC.gp.gp_Pnt(border_pnt[0], border_pnt[1], max_z)
            # Create wire (border of one face)
            wire = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeWire()
            # 0
            edge0 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point0, point1)
            wire.Add(edge0.Edge())
            cap1_wire.Add(edge0.Edge())
            # 1
            edge1 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point1, point2)
            wire.Add(edge1.Edge())
            # 2
            edge2 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point2, point3)
            wire.Add(edge2.Edge())
            cap2_wire.Add(edge2.Edge())
            # 3
            edge3 = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeEdge(point3, point0)
            wire.Add(edge3.Edge())
            # Create face from the wire
            face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(wire.Wire())
            sewing.Add(face.Shape())

            cap1_face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(cap1_wire.Wire())
            sewing.Add(cap1_face.Shape())
            cap2_face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(cap2_wire.Wire())
            sewing.Add(cap2_face.Shape())

            # Sew it all together
            sewing.Perform()
            sewing_shape = sewing.SewedShape()
            # Create shell, solid and compound
            shell = OCC.TopoDS.topods_Shell(sewing_shape)
            make_solid = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeSolid()
            make_solid.Add(shell)
            solid = make_solid.Solid()
            builder = OCC.BRep.BRep_Builder()
            builder.MakeSolid(solid)
            builder.Add(solid, shell)

            # Try to fix orientation of solid volume
            OCC.BRepLib.breplib().OrientClosedSolid(solid)

            self.area_volumes[border_id] = solid

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

            # Use OCC to approximate rivers with BSpline curve
            array = OCC.TColgp.TColgp_Array1OfPnt(1, len(self.rivers_data_3d[river_id]))
            point_id = 1
            for point in self.rivers_data_3d[river_id]:
                array.SetValue(point_id, OCC.gp.gp_Pnt(point[0], point[1], point[2]))
                point_id += 1
            self.rivers_curves_3d[river_id] = OCC.GeomAPI.GeomAPI_PointsToBSpline(array).Curve()

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

        diffs = self.conf['terrain']['approximation']['differences']
        display_surf = self.conf['display']['surface']
        display_terr = self.conf['display']['terrain']
        display_rivers = self.conf['display']['rivers']
        display_fractures = self.conf['display']['fractures']

        if display_surf is True:
            if self.volume is not None:
                # Display IDs of points at are border
                self.approximate_2d_borders()
                # for point_id, point in enumerate(self.area_borders_3d[0]):
                #     occ_point = OCC.gp.gp_Pnt(point[0], point[1], point[2] + 10.0)
                #     display.DisplayMessage(occ_point, str(point_id))
                # Display volume
                display.DisplayShape(self.volume, update=True)
                # Display original terrain with transparency
                if self.shell is not None:
                    ais_shell = display.DisplayShape(self.shell, update=True)
                    display.Context.SetTransparency(ais_shell, 0.8)
                else:
                    ais_sewing = display.DisplayShape(self.sewing_shape, update=True)
                    display.Context.SetTransparency(ais_sewing, 0.8)
            elif self.shell is not None:
                display.DisplayShape(self.shell, update=True)
            else:
                # Draw all bspline surfaces
                for bspline in self.bspline_surfaces.values():
                    display.DisplayShape(bspline.GetHandle(), update=True)

        if display_terr is True:
            # Draw points of terrain
            a_presentation, group = self.create_ogl_group(display)
            black = OCC.Quantity.Quantity_Color(OCC.Quantity.Quantity_NOC_BLACK)
            asp = OCC.Graphic3d.Graphic3d_AspectLine3d(black, OCC.Aspect.Aspect_TOL_SOLID, 1)

            if diffs is True:
                pnt_array = OCC.Graphic3d.Graphic3d_ArrayOfPoints(self.point_count,
                                                                  True,  # hasVColors
                                                                  )
            else:
                pnt_array = OCC.Graphic3d.Graphic3d_ArrayOfPoints(self.point_count,
                                                                  False,  # hasVColors
                                                                  )
            if diffs is True:
                idx = 1
            # Default RGB color of point (white)
            for point in self.terrain_data:
                pnt = OCC.gp.gp_Pnt(point[0], point[1], point[2])
                pnt_array.AddVertex(pnt)
                if diffs is True:
                    # create the point, with a diff color
                    if self.max_diff > 0.0:
                        diff = self.tW[idx - 1] / self.max_diff
                    else:
                        diff = 0.0
                    rgb = colorsys.hsv_to_rgb(diff, 1.0, 1.0)
                    pnt_array.SetVertexColor(idx, rgb[0], rgb[1], rgb[2])
                    idx += 1

            group.SetPrimitivesAspect(asp.GetHandle())
            group.AddPrimitiveArray(pnt_array.GetHandle())
            a_presentation.Display()

        if display_rivers is True:
            a_presentation, group = self.create_ogl_group(display)
            black = OCC.Quantity.Quantity_Color(OCC.Quantity.Quantity_NOC_BLACK)
            asp = OCC.Graphic3d.Graphic3d_AspectLine3d(black, OCC.Aspect.Aspect_TOL_SOLID, 1)
            count = 0
            for river in self.rivers_data_3d.values():
                count += len(river)
            pnt_array = OCC.Graphic3d.Graphic3d_ArrayOfPoints(count,
                                                              False,  # hasVColors
                                                              )
            for river in self.rivers_data_3d.values():
                for point in river:
                    pnt = OCC.gp.gp_Pnt(point[0], point[1], point[2])
                    pnt_array.AddVertex(pnt)
            group.SetPrimitivesAspect(asp.GetHandle())
            group.AddPrimitiveArray(pnt_array.GetHandle())
            a_presentation.Display()

            # for river in self.rivers_curves_3d.values():
            #     display.DisplayShape(river, update=True)

        if display_fractures is True and 'input' in self.conf['fractures']:
            display.DisplayShape(self.fractures_shape, update=True)

        start_display()

    @staticmethod
    def cap_extrude(bspl_surf, z_diff):
        """
        Cap extrusion
        :param bspl_surf:
        :param z_diff:
        :return OCC B_Spline
        """
        # Non-periodic surface
        uperiod = False
        vperiod = False

        udeg = bspl_surf.UDegree()
        vdeg = bspl_surf.VDegree()

        nb_u_poles = bspl_surf.NbUPoles()
        nb_v_poles = bspl_surf.NbVPoles()

        poles = OCC.TColgp.TColgp_Array2OfPnt(1, nb_u_poles, 1, nb_v_poles)
        bspl_surf.Poles(poles)
        new_poles = OCC.TColgp.TColgp_Array2OfPnt(1, nb_u_poles, 1, nb_v_poles)

        nb_u_knots = nb_u_mults = bspl_surf.NbUKnots()
        nb_v_knots = nb_v_mults = bspl_surf.NbVKnots()
        uknots = OCC.TColStd.TColStd_Array1OfReal(1, nb_u_knots)
        vknots = OCC.TColStd.TColStd_Array1OfReal(1, nb_v_knots)
        umults = OCC.TColStd.TColStd_Array1OfInteger(1, nb_u_mults)
        vmults = OCC.TColStd.TColStd_Array1OfInteger(1, nb_v_mults)
        bspl_surf.UKnots(uknots)
        bspl_surf.VKnots(vknots)
        bspl_surf.VMultiplicities(umults)
        bspl_surf.VMultiplicities(vmults)

        # Set Z coordinate in poles
        for u_idx in range(1, nb_u_poles + 1):
            for v_idx in range(1, nb_v_poles + 1):
                _u_idx = nb_u_poles + 1 - u_idx
                _v_idx = nb_v_poles + 1 - v_idx
                point = poles.Value(_u_idx, _v_idx)
                new_poles.SetValue(u_idx, v_idx, OCC.gp.gp_Pnt(point.X(), point.Y(), point.Z() + z_diff))

        return OCC.Geom.Geom_BSplineSurface(new_poles, uknots, vknots, umults, vmults, udeg, vdeg, uperiod, vperiod)

    @staticmethod
    def extrude_edge(upole, uknots, umult, udeg, z_diff):
        """
        Extrude one edge.
        :param upole:
        :param uknots:
        :param umult:
        :param udeg:
        :param z_diff:
        :return OCC B-Spline
        """

        # Non-periodic surface
        uperiod = False
        vperiod = False

        vdeg = 1

        # Create 2D array of poles (control points)
        poles = OCC.TColgp.TColgp_Array2OfPnt(1, upole.Length(), 1, 2)
        for index in range(1, upole.Length() + 1):
            point = upole.Value(index)
            poles.SetValue(index, 1, point)
            poles.SetValue(index, 2, OCC.gp.gp_Pnt(point.X(), point.Y(), point.Z() + z_diff))

        # Length of uknots and umult has to be same
        # Same rule is for vknots and vmult
        vknot_len = vmult_len = 2

        # Knots for V direction
        vknots = OCC.TColStd.TColStd_Array1OfReal(1, vknot_len)

        # Main curves begins and ends at first and last points
        vknots.SetValue(1, 0.0)
        vknots.SetValue(2, 1.0)

        # Multiplicities for U and V direction
        vmult = OCC.TColStd.TColStd_Array1OfInteger(1, vmult_len)

        # First and last multiplicities are set to udeg + 1 (vdeg respectively),
        # because we want main curves to start and finish on the first and
        # the last points
        vmult.SetValue(1, vdeg + 1)
        vmult.SetValue(2, vdeg + 1)

        # Try to create surface
        return OCC.Geom.Geom_BSplineSurface(poles, uknots, vknots, umult, vmult, udeg, vdeg, uperiod, vperiod)

    def extrude_surface(self, bspl_surf, z_diff=-100):
        """
        This method extrude OCC B-Spline surface and it returns list of B-Spline surfaces. These surfaces can be
        assembled into the volume object.
        :param bspl_surf:
        :param z_diff:
        :return: list of B-Spline surfaces
        """

        surfaces = []

        nb_u_poles = bspl_surf.NbUPoles()
        nb_v_poles = bspl_surf.NbVPoles()

        poles = OCC.TColgp.TColgp_Array2OfPnt(1, nb_u_poles, 1, nb_v_poles)
        bspl_surf.Poles(poles)

        nb_u_knots = nb_u_mults = bspl_surf.NbUKnots()
        uknots = OCC.TColStd.TColStd_Array1OfReal(1, nb_u_knots)
        umults = OCC.TColStd.TColStd_Array1OfInteger(1, nb_u_mults)
        bspl_surf.UKnots(uknots)
        bspl_surf.UMultiplicities(umults)

        # 1: Extrude one border edge
        edge = OCC.TColgp.TColgp_Array1OfPnt(1, nb_u_poles)
        for index in range(1, nb_u_poles + 1):
            _index = nb_u_poles + 1 - index
            edge.SetValue(index, poles.Value(1, _index))
        extruded_surf = self.extrude_edge(edge, uknots, umults, bspl_surf.UDegree(), z_diff)
        surfaces.append(extruded_surf)

        # 2: Extrude one border edge
        edge = OCC.TColgp.TColgp_Array1OfPnt(1, nb_u_poles)
        for index in range(1, nb_u_poles + 1):
            edge.SetValue(index, poles.Value(nb_v_poles, index))
        extruded_surf = self.extrude_edge(edge, uknots, umults, bspl_surf.UDegree(), z_diff)
        surfaces.append(extruded_surf)

        # 3: Extrude one border edge
        edge = OCC.TColgp.TColgp_Array1OfPnt(1, nb_v_poles)
        for index in range(1, nb_v_poles + 1):
            edge.SetValue(index, poles.Value(index, 1))
        extruded_surf = self.extrude_edge(edge, uknots, umults, bspl_surf.UDegree(), z_diff)
        surfaces.append(extruded_surf)

        # 4: Extrude one border edge
        edge = OCC.TColgp.TColgp_Array1OfPnt(1, nb_v_poles)
        for index in range(1, nb_v_poles + 1):
            _index = nb_v_poles + 1 - index
            edge.SetValue(index, poles.Value(_index, nb_u_poles))
        extruded_surf = self.extrude_edge(edge, uknots, umults, bspl_surf.UDegree(), z_diff)
        surfaces.append(extruded_surf)

        # Cap it
        cap_surf = self.cap_extrude(bspl_surf, z_diff)
        surfaces.append(cap_surf)

        return surfaces

    def output_approx_data(self):
        """
        Try to output approximated data to BREP file format
        """

        extrude_diff = self.conf['terrain']['extrude_diff']

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

        if extrude_diff != 0.0:
            # TODO: this code works properly only in situation, when surface is approximated only with one B-Spline
            surfaces = self.extrude_surface(occ_bspline, extrude_diff)
            # Added extruded surfaces to the sewing
            for surface in surfaces:
                face = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeFace(surface.GetHandle(), error).Shape()
                sewing.Add(face)
            # Sew it all together
            sewing.Perform()
            self.sewing_shape = sewing.SewedShape()
            # Create shell, solid and compound
            self.shell = OCC.TopoDS.topods_Shell(self.sewing_shape)
            make_solid = OCC.BRepBuilderAPI.BRepBuilderAPI_MakeSolid()
            make_solid.Add(self.shell)
            solid = make_solid.Solid()
            # Try to fix solid shape
            OCC.BRepLib.breplib().OrientClosedSolid(solid)
            builder = OCC.BRep.BRep_Builder()
            builder.MakeSolid(solid)
            builder.Add(solid, self.shell)
            compound = OCC.TopoDS.TopoDS_Compound()
            builder = OCC.BRep.BRep_Builder()
            builder.MakeCompound(compound)
            if len(self.area_borders_2d) > 0:
                self.create_volume_from_area(extrude_diff)
                # TODO: This code works only for one area
                print('Computing union between volumes ...')
                start_time = time.time()
                self.volume = OCC.BRepAlgoAPI.BRepAlgoAPI_Common(self.area_volumes[0], solid).Shape()
                end_time = time.time()
                print('Computed in {0} seconds.'.format(end_time - start_time))
                builder.Add(compound, self.volume)
            else:
                builder.Add(compound, solid)
            # Write compound to the BREP file
            OCC.BRepTools.breptools_Write(compound, self.conf['output'])
        else:
            sewing.Perform()
            self.sewing_shape = sewing.SewedShape()
            if len(self.area_borders_2d) > 0:
                self.create_volume_from_area(extrude_diff)
                # Test of section between surface and faces
                # print('Computing section between volume and surface ...')
                # start_time = time.time()
                # OCC.BRepAlgoAPI.BRepAlgoAPI_Section(self.area_volumes[0], self.sewing_shape)
                # end_time = time.time()
                # print('Computed in {0} seconds.'.format(end_time - start_time))
                # TODO: This code works only for one area
                print('Computing union between volume and surface ...')
                start_time = time.time()
                self.volume = OCC.BRepAlgoAPI.BRepAlgoAPI_Common(self.area_volumes[0], self.sewing_shape).Shape()
                end_time = time.time()
                print('Computed in {0} seconds.'.format(end_time - start_time))
                OCC.BRepTools.breptools_Write(self.volume, self.conf['output'])
            else:
                OCC.BRepTools.breptools_Write(self.sewing_shape, self.conf['output'])
