"""
This module is used for loading terrain data
"""

import yaml
import math
import sys

from scipy import interpolate
import numpy as np

import convert.bspline

class TerrainData(object):
    """
    Class representing Terrain data (terrain, rivers, aproximation, etc)
    """
    def __init__(self, yaml_file_name):
        """Constructor of TerrainData"""
        super(TerrainData, self).__init__()
        self.yaml_file_name = yaml_file_name
        self.conf = {} # Configuration loaded from yaml file
        self.terrain_data = [] # Terrain data list of tuples (3d coordinates)
        self.point_count = 0
        self.tck = {} # Patches of b-spline surface
        self.tX = None #
        self.tY = None # Terrain X,Y,Z cache for numpy API
        self.tZ = None #
        self.tW = None # Visualization purpose
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
        self.grid = {} # Cache for bilinear interpolation
        self.rivers_data_2d = {}
        self.rivers_data_3d = {}
        self.area_borders_2d = {}
        self.area_borders_3d = {}

    def load_conf_from_yaml(self):
        """
        Load configuration form yaml file
        """
        with open(self.yaml_file_name, 'r') as yaml_file:
            self.conf = yaml.load(yaml_file)

    def __post_process_terrain_data(self):
        """
        Try to postprocess terrain data
        """
        self.point_count = len(self.terrain_data)
        self.tX = [item[0] for item in self.terrain_data]
        self.tY = [item[1] for item in self.terrain_data]
        self.tZ = [item[2] for item in self.terrain_data]

        def find_first_different(array):
            """Find index of first different item"""
            index = 1
            first = array[0]
            while array[index] == first:
                index += 1
            index -= 1
            if index < len(array):
                return index
            else:
                return -1

        def find_first_same(array):
            """Find index of first same item"""
            index = 1
            first = array[0]
            while array[index] != first:
                index += 1
            index -= 1
            if index < len(array):
                return index
            else:
                return -1

        self.min_x = min(self.tX)
        self.max_x = max(self.tX)
        self.min_y = min(self.tY)
        self.max_y = max(self.tY)
        # Try to compute size of 2d array accordin repeated values
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
        for index in range(0, self.point_count):
            i = int( math.floor( (self.terrain_data[index][0] - self.min_x) / self.dx ) )
            j = int( math.floor( (self.terrain_data[index][1] - self.min_y) / self.dy ) )
            self.grid[(i, j)] = (self.terrain_data[index][0], self.terrain_data[index][1], self.terrain_data[index][2])
        self.bspline_surfaces = {}

    def load_terrain(self):
        """
        Try to load data of terain
        """
        with open(self.conf['terrain'], 'r') as data_file:
            for line in data_file:
                self.terrain_data.append( tuple( float(item) for item in line.split() ) )
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
                self.area_borders_2d[0].append( tuple( float(item) for item in items[1:]) )

    def aproximate_terrain(self):
        """
        Try to aproximate terrain with bspline surface
        """
        tck,fp,ior,msg = interpolate.bisplrep(self.tX, self.tY, self.tZ, kx=5, ky=5, full_output=1)
        self.tck[(self.min_x, self.min_y, self.max_x, self.max_y)] = tck
        # Compute difference between original terrain data and b-spline surface
        self.tW = [abs(it[2] - interpolate.bisplev(it[0], it[1], tck)) for it in self.terrain_data]

    def aproximate_2d_borders(self):
        """
        Try to aproximate z coordinates of borders using terrain data
        """
        
        for border_id,border_2d in self.area_borders_2d.items():
            self.area_borders_3d[border_id] = []
            for bp in border_2d:
                # Compute indexes to the grid first
                i = int(math.floor( (bp[0] - self.min_x) / self.dx ))
                j = int(math.floor( (bp[1] - self.min_y) / self.dy ))
                # Compute weights for bilineral intebpolation
                kx = (bp[0] - (self.min_x + self.dx * i)) / self.dx
                ky = (bp[1] - (self.min_y + self.dy * j)) / self.dy
                z1 = self.grid[(i, j)][2]
                z2 = self.grid[(i + 1, j)][2]
                z3 = self.grid[(i, j + 1)][2]
                z4 = self.grid[(i + 1, j + 1)][2]
                z12 = (1.0 - kx) * z1 + kx * z2
                z34 = (1.0 - kx) * z3 + kx * z4
                Z = (1.0 - ky) * z12 + ky * z34
                self.area_borders_3d[border_id].append( (bp[0], bp[1], Z) )


    def aproximate_2d_rivers(self):
        """
        Try to aproximate z coordinates of rivers using terrain data
        """
        
        for river_id,river2d in self.rivers_data_2d.items():
            self.rivers_data_3d[river_id] = []
            last_z = sys.maxsize
            for rp in river2d:
                # Compute indexes to the grid first
                i = int(math.floor( (rp[0] - self.min_x) / self.dx ))
                j = int(math.floor( (rp[1] - self.min_y) / self.dy ))
                # Compute weights for bilineral interpolation
                kx = (rp[0] - (self.min_x + self.dx * i)) / self.dx
                ky = (rp[1] - (self.min_y + self.dy * j)) / self.dy
                z1 = self.grid[(i, j)][2]
                z2 = self.grid[(i + 1, j)][2]
                z3 = self.grid[(i, j + 1)][2]
                z4 = self.grid[(i + 1, j + 1)][2]
                z12 = (1.0 - kx) * z1 + kx * z2
                z34 = (1.0 - kx) * z3 + kx * z4
                Z = (1.0 - ky) * z12 + ky * z34
                # Cut out too big z values
                if last_z < Z:
                    #print 'last_z: ', last_z, ' < Z: ', Z, ' dz: ', Z - last_z
                    Z = last_z
                last_z = Z 
                self.rivers_data_3d[river_id].append( (rp[0], rp[1], Z) )

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
        start_display()

    def output_approx_data(self):
        """
        Try to output approximated data to BREP file format
        """

        for key,scipy_bspline in self.tck.items():
            occ_bspline = convert.bspline.scipy_to_occ(scipy_bspline)
            self.bspline_surfaces[key] = occ_bspline