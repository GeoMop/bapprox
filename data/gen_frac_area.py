"""
Try to generate fractal area
"""

import argparse
import random
import matplotlib.pyplot as plt


def load_area(filename):
    """
    Try to load definition of area border
    """
    border = []
    with open(filename, 'r') as area_file:
        for line in area_file:
            items = line.split()
            # Only one border ATM
            border.append(list(float(item) for item in items[1:]))
    return border


def write_area(border, filename):
    """
    Try to write border to filename
    :param border:
    :param filename:
    :return:
    """
    with open(filename, 'w') as output_file:
        for point_id, point in enumerate(border):
            output_file.write("{0} {1} {2} {3}\n".format(
                point_id,
                point[0],
                point[1],
                point[2]
            ))


def new_point(prev_point, point, factor):
    """
    Try to generate new point from two points
    :param prev_point:
    :param point:
    :param factor:
    :return:
    """
    mid_point = [0.0, 0.0, 0.0]
    # X
    dx = point[0] - prev_point[0]
    diff_x = factor * (random.random() - 0.5)
    mid_point[0] = prev_point[0] + dx / 2.0 + diff_x
    # Y
    dy = point[1] - prev_point[1]
    diff_y = factor * (random.random() - 0.5)
    mid_point[1] = prev_point[1] + dy / 2.0 + diff_y
    # Z
    mid_point[2] = prev_point[2] + (point[2] - prev_point[2]) / 2.0
    return mid_point


def frac_border(border, factor, level=1):
    """
    Frac border
    :param border:
    :param factor:
    :param level
    :return:
    """

    points = []
    prev_point = None
    first_point = None
    for point in border:
        if prev_point is not None:
            mid_point = new_point(prev_point, point, factor)
            points.append(mid_point)
        else:
            first_point = point
        points.append(point)
        prev_point = point
    mid_point = new_point(point, first_point, factor)
    points.append(mid_point)
    if level == 0:
        return points
    else:
        return frac_border(points, factor/2.0, level-1)


def draw_border(border, color='k', show=True):
    """
    Try to draw border as polyline
    :param border:
    :param color:
    :param show:
    :return:
    """
    prev_point = None
    first_point = None

    x_coords = [point[0] for point in border]
    y_coords = [point[1] for point in border]

    plt.plot(x_coords, y_coords, "o", color=color)

    for point in border:
        if prev_point is not None:
            x_coords = [prev_point[0], point[0]]
            y_coords = [prev_point[1], point[1]]
            plt.plot(x_coords, y_coords, color=color)
        else:
            first_point = point
        prev_point = point
    x_coords = [first_point[0], point[0]]
    y_coords = [first_point[1], point[1]]
    plt.plot(x_coords, y_coords, color=color)
    if show is True:
        plt.show()


def solve(input_filename, output_filename, factor, level):
    """
    My function of this script
    :param input_filename:
    :param output_filename:
    :param factor:
    :param level:
    :return:
    """
    border = load_area(input_filename)
    draw_border(border, color='r', show=False)
    new_border = frac_border(border, factor, level)
    draw_border(new_border)
    write_area(new_border, output_filename)


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="Yaml conf file", default=None)
    parser.add_argument("-o", "--output", type=str, required=False,
                        help="Output file", default="output.txt")
    parser.add_argument("-l", "--level", type=int, required=False,
                        help="Level of recursive", default=0)
    parser.add_argument("-f", "--factor", type=int, required=False,
                        help="Factor of change", default=100)
    args = parser.parse_args()
    return args.input, args.output, args.factor, args.level


def main():
    """
    Main function
    :return:
    """
    input_filename, output_filename, factor, level = parse_arguments()
    solve(input_filename, output_filename, factor, level)


if __name__ == '__main__':
    main()
