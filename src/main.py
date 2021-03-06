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

    if yaml_conf_filename is None:
        return

    terrain = terrain_data.TerrainData(yaml_conf_filename)

    terrain.load_conf_from_yaml()
    terrain.load_terrain()
    terrain.load_rivers()
    terrain.load_area()
    terrain.load_fractures()

    terrain.approximate_terrain()

    terrain.approximate_2d_rivers()

    terrain.output_approx_data()

    if terrain.conf['display']['terrain'] is True or terrain.conf['display']['surface'] is True:
        terrain.display_terrain()


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filename", type=str, required=True,
                        help="Yaml conf file", default=None)
    args = parser.parse_args()
    return args.filename

if __name__ == '__main__':
    # Parse argument
    main(parse_arguments())
