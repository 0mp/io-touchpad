#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""io-touchpad application.

This application will learn to associate symbols you draw on your touchpad
with certain actions you bound to those symbols.

Your finger's absolute position is fetched using a library called touchpadlib
which is based on evtest (version 1.32).

Source code encoding: UTF-8
"""

import queue
import sys

from terminationhandler import terminationhandler
from threads import application
from threads import listener
from touchpadlib import touchpadlib
from classifier import classifier

def main():
    """The main function."""
    # SIGINT signal handler.
    terminationhandler.setup()

    training_size = 0
    learning_mode = len(sys.argv) > 1
    if learning_mode:
        try:
            training_size = int(sys.argv[1])
        except ValueError:
            print("wrong argument")
            sys.exit(1)
        if training_size == 0:
            clsf = classifier.Classifier()
            clsf.learn()
            sys.exit(0)
        elif training_size < 3:
            print("argument is to little")
            sys.exit(1)

    thread_queue = queue.Queue()

    touchpad_specification = touchpadlib.Touchpadlib.get_specification()

    print("Touchpad specification: "
          "min/max x: %d/%d, min/max y: %d/%d, min/max pressure %d/%d." %
          (touchpad_specification['min_x'],
           touchpad_specification['max_x'],
           touchpad_specification['min_y'],
           touchpad_specification['max_y'],
           touchpad_specification['min_pressure'],
           touchpad_specification['max_pressure']))

    # Greating.
    if learning_mode:
        print("Welcome in learning mode. Please draw %d versions of the symbol." % (training_size))
    else:
        print("Use your touchpad as usual. Have a nice day!")

    # Run both threads.
    listener.start(thread_queue)
    application.application_thread(thread_queue, learning_mode, training_size)

main()
