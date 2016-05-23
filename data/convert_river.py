"""
Try to generate fractal area
"""

import argparse
import math
import matplotlib.pyplot as plt


def load_rivers(filename):
    """
    Try to load definition of area border
    :param filename: Filename of intput file
    :return: Dictionary of list containing list of river points
    """
    rivers = {}
    with open(filename, 'r') as area_file:
        for line in area_file:
            items = line.split()
            river_id = int(items[3])
            try:
                river = rivers[river_id]
            except KeyError:
                river = rivers[river_id] = []
            # Only one border ATM
            river.append([float(items[0]), float(items[1]), float(items[2])])
    return rivers


def write_rivers(rivers, filename):
    """
    Try to write rivers to filename
    :param rivers: Dictionary of rivers, each river is list of river coordinate
    :param filename: Filename of output file
    :return: None
    """
    with open(filename, 'w') as output_file:
        for river_id, river in rivers.items():
            for point in river:
                output_file.write("{0} {1} {2} {3}\n".format(
                    point[0],
                    point[1],
                    point[2],
                    river_id
                ))


def draw_rivers(rivers, color='k'):
    """
    Try to draw river as polyline
    :param rivers: Dictionary of rivers, each river is list of river coordinate
    :param color: Color of rivers
    :return: None
    """

    for river in rivers.values():
        prev_point = None
        x_coords = [point[0] for point in river]
        y_coords = [point[1] for point in river]

        plt.plot(x_coords, y_coords, "o", color=color)

        for point in river:
            if prev_point is not None:
                x_coords = [prev_point[0], point[0]]
                y_coords = [prev_point[1], point[1]]
                plt.plot(x_coords, y_coords, color=color)
            prev_point = point
        plt.plot(x_coords, y_coords, color=color)
    plt.show()


def merge_points(river, max_diff):
    """
    This function tries to merge two points too near together.
    :param river: List of river points.
    :param max_diff: Maximal distance between two points of river
    :return: Modified list of river points.
    """
    new_river = []
    prev_point = None
    for point in river:
        if prev_point is not None:
            diff_x = (point[0] - prev_point[0])
            diff_y = (point[1] - prev_point[1])
            diff = math.sqrt(diff_x**2 + diff_y**2)
            if diff < max_diff:
                # Ignore point
                pass
            else:
                # Add point to the new list
                new_river.append(point)
        prev_point = point
    return new_river


def split_river(river, max_diff):
    """
    Try to split points of river to be equidistant
    :param river: List of river points.
    :param max_diff: Maximal distance between two points of river
    :return: Modified list of river points.
    """
    new_river = []
    prev_point = None
    for point in river:
        if prev_point is not None:
            diff_x = (point[0] - prev_point[0])
            diff_y = (point[1] - prev_point[1])
            diff = math.sqrt(diff_x**2 + diff_y**2)
            if diff > max_diff:
                total_dist_01 = point[2] - prev_point[2]
                num_cuts = diff / max_diff
                fraction = math.modf(num_cuts)
                if fraction < 0.5:
                    num_cuts = math.floor(num_cuts)
                else:
                    num_cuts = math.ceil(num_cuts)
                diff_x = (point[0] - prev_point[0]) / num_cuts
                diff_y = (point[1] - prev_point[1]) / num_cuts
                diff_dist_01 = (point[2] - prev_point[2]) / num_cuts
                coord_x = prev_point[0] + diff_x
                coord_y = prev_point[1] + diff_y
                dist_01 = diff_dist_01
                while dist_01 < total_dist_01:
                    new_river.append([coord_x, coord_y, prev_point[2] + dist_01])
                    coord_x += diff_x
                    coord_y += diff_y
                    dist_01 += diff_dist_01
        new_river.append(point)
        prev_point = point
    return new_river


def convert(input_filename, output_filename, limit, diff):
    """
    Function converting river points to be more than less equidistant
    :param input_filename: Filename of input file
    :param output_filename: Filename of output file
    :param limit: Minimal distance between two points
    :param diff: Maximal distance between two points
    :return: None
    """
    rivers = load_rivers(input_filename)
    # TODO: convert rivers
    new_rivers = {}
    for river_id, river in rivers.items():
        river = split_river(river, diff)
        if limit > 0.0:
            river = merge_points(river, limit)
        new_rivers[river_id] = river
    draw_rivers(new_rivers, color='r')
    write_rivers(new_rivers, output_filename)


def parse_arguments():
    """
    Parse command line arguments
    :return: tuple of values: (input filename, output filename, limit, )
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Yaml conf file", default=None)
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Output file", default="output.txt")
    parser.add_argument("-l", "--limit", type=float, required=False,
                        help="Limit maximal distance between two points (near points are merged).", default=0.0)
    parser.add_argument("-d", "--diff", type=float, required=False,
                        help="Minimal difference between two points (too far points are divided to several points)",
                        default=1.0)
    args = parser.parse_args()
    return args.input, args.output, args.limit, args.diff


def main():
    """
    Main function
    :return: None
    """
    input_filename, output_filename, limit, diff = parse_arguments()
    convert(input_filename, output_filename, limit, diff)


if __name__ == '__main__':
    main()
