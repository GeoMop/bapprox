#!/usr/bin/env python

"""
This application can be used for approximation of surface using
bspline surfaces.
"""

import argparse
import terrain_data


def main(yaml_conf_filename):
    """
    Try to read data from all files
    :param yaml_conf_filename: dictionary with configuration
    """
    terrain = terrain_data.TerrainData(yaml_conf_filename)

    terrain.load_conf_from_yaml()
    terrain.load_terrain()
    terrain.load_rivers()
    terrain.load_area()

    terrain.approximate_terrain()

    if terrain.conf['area']['approximate'] is True:
        terrain.approximate_2d_borders()

    if terrain.conf['rivers']['approximate'] is True:
        terrain.approximate_2d_rivers()

    terrain.output_approx_data()

    if terrain.conf['display']['terrain'] is True or terrain.conf['display']['surface']:
        terrain.display_terrain()


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", type=str,
                        help="Yaml conf file", default=None)
    args = parser.parse_args()
    main(args.filename)

if __name__ == '__main__':
    # Parse argument
    parse_arguments()
