#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""This is a utility to visualize the symbols drawn on touchpads.

It starts its own modified version of the app/app.py and creates a figure
for every signal_list the modified thread/application.py would send to the
interpreter.

The touchpad size is detected automatically. A suitable tolerance might be
added if the detected size is leading to the out of bounds errors.

You'll need matplotlib to make it work.
Sometimes matplotlib doesn't work when using virtual environment
like pyvenv. It is advised to download it using your system
package manager.

@TODO I don't know how to deal with the relative paths in Python3.
      This is why we have to run the program from the app/ directory and
      use hacks like `sys.path.append(".")` to this program work. (~Mateusz)
"""

import argparse
import queue
import time
import types
import sys
import matplotlib.pyplot as plt

sys.path.append(".")

from threads import listener, application
from terminationhandler import terminationhandler
from touchpadlib import touchpadlib


class MatrixAnalyser(object):
    """The matrix analyser class.

    The main purpose of this class is to store the information about the
    touchpad and be able to use them inside save_symbols() after the
    application.send_points_to_interpreter() function is monkey-patch with the
    save_symbols() function.

    Attributes:
        min_x (int): The minimal value of the x axis on the touchpad.
        max_x (int): The maximal value of the x axis on the touchpad.
        min_y (int): The minimal value of the y axis on the touchpad.
        max_y (int): The maximal value of the y axis on the touchpad.
        tolerance(int): The parameter added to the maximal values to increase
            the range of the touchpad.

    Constants:
        FIGURES_PATH (str): The path to the directory where the figures are
            saved.
        FIGURE_DPI (int): The quality of the output figures.
    """

    FIGURES_PATH = "./tools/data/matrixanalyser/figures/"
    FIGURE_DPI = 500

    def __init__(self, specification, tolerance):
        """Constructor.

        Args:
            specification (dict): The dictionary with the specification of
                the touchpad.
            tolerance(int): The parameter added to the maximal values to
                increase the range of the touchpad.
        """
        self.min_x = specification['min_x']
        self.max_x = specification['max_x'] + tolerance
        self.min_y = specification['min_y']
        self.max_y = specification['max_y'] + tolerance
        self.tolerance = tolerance

    def save_symbols(self, signal_list):
        """Generate figures based on the signal_list.

        This is a function which monkey-patches the
        application.send_points_to_interpreter() function. It saves generated
        figures to the FIGURES_PATH. It doesn't save the figure if one of the
        signal went out of bound or if the signal_list is empty.

        Args:
            signal_list (list): List of TouchpadSingals received from the
                touchpad.
        """
        if not signal_list:
            return

        x_coordinates = []
        y_coordinates = []

        for signal in signal_list:
            normalized_x = self._normalize_x(signal.get_x())
            if normalized_x is None:
                return
            x_coordinates.append(normalized_x)

            normalized_y = self._normalize_y(signal.get_y())
            if normalized_y is None:
                return
            y_coordinates.append(normalized_y)

        plt.gca().invert_yaxis()
        plt.scatter(x_coordinates, y_coordinates, s=3)
        name = "figure-" + str(int(time.time())) + '.png'
        plt.savefig(self.FIGURES_PATH + name, dpi=self.FIGURE_DPI)
        #  plt.show()

    def _normalize_x(self, x_value):
        normalized_x = x_value - self.min_x
        if normalized_x > self.max_x:
            print("matrixanalyser.py: warning: "
                  "The x coordinate of the recived signal "
                  "is out of the touchpad's detected range. "
                  "Consider passing a larger tolerance argument.",
                  file=sys.stderr)
            return None
        else:
            return normalized_x

    def _normalize_y(self, y_value):
        normalized_y = y_value - self.min_y
        if normalized_y > self.max_y:
            print("matrixanalyser.py: warning:"
                  "The y coordinate of the recived signal "
                  "is out of the touchpad's detected range. "
                  "Consider passing a larger tolerance argument.",
                  file=sys.stderr)
            return None
        else:
            return normalized_y


def main():
    """The main function."""
    terminationhandler.setup()

    description = 'Visualize symbols.'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--tolerance', dest='tolerance', default=0,
                        help="set a tolerance for the resolution "
                        "(default: 0)")
    args = parser.parse_args()

    touchpad_specification = touchpadlib.Touchpadlib.get_specification()

    thread_queue = queue.Queue()

    matrix_analyser = MatrixAnalyser(touchpad_specification, args.tolerance)
    application.send_points_to_interpreter = \
        types.MethodType(MatrixAnalyser.save_symbols, matrix_analyser)

    listener.start(thread_queue)
    application.application_thread(thread_queue)

main()
