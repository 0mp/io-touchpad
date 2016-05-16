#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""io-touchpad application.

This application will learn to associate symbols you draw on your touchpad
with certain actions you bound to those symbols.

Your finger's absolute position is fetched using a library called touchpadlib
which is based on evtest (version 1.32).

Source code encoding: UTF-8
"""

import argparse
import queue
import sys

from terminationhandler import terminationhandler
from threads import application
from threads import listener
from classifier import classifier

MIN_TRAINING_SIZE = 3


def _get_configured_parser():
    """Configure the commandline arguments parser."""
    description = 'Teaching your touchpad magic tricks since 2016.'
    parser = argparse.ArgumentParser(add_help=True, description=description)

    group = parser.add_mutually_exclusive_group()

    group.add_argument('-b', '--bitness', dest='system_bitness',
                       default=None, metavar='BITS',
                       help='set the bitness of your operation system; '
                       'this option triggers the use of the hardcoded '
                       'symbols', choices={'32', '64'})
    group.add_argument('-l', '--learning', dest='training_size',
                       default=None, metavar='SIZE',
                       help='start the application in the learning mode; '
                       'you will be asked to draw a symbol SIZE times; '
                       'SIZE should be at least ' + str(MIN_TRAINING_SIZE),
                       type=int)
    group.add_argument('-r', '--repeat', dest='repeat_classification',
                       default=False, action='store_true',
                       help='repeat the classification on the latest '
                       'user-defined set of drawings')

    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument('-s', '--symbol', dest='symbol_name',
                       default=None, metavar='NAME',
                       help='name of symbol to learn in the learning mode')
    return parser


def main():
    """The main function."""
    # SIGINT signal handler.
    terminationhandler.setup()

    parser = _get_configured_parser()
    args = parser.parse_args()

    training_size = args.training_size
    learning_mode = training_size is not None
    symbol_name = args.symbol_name

    if args.repeat_classification:
        print('Repeating the classification within the learning process.')
        clsf = classifier.Classifier(True)
        clsf.learn(True, symbol_name)
        sys.exit(0)

    if learning_mode:
        print(args.symbol_name)
        if training_size < 3:
            print('app.py: error: the training size should be at least %d'
                  % (MIN_TRAINING_SIZE), file=sys.stderr)
            sys.exit(1)

    system_bitness = args.system_bitness
    if system_bitness is not None:
        system_bitness = int(system_bitness)

    thread_queue = queue.Queue()

    # Run both threads.
    listener.start(thread_queue)
    application.application_thread(thread_queue, learning_mode, training_size,
                                   system_bitness, symbol_name)

main()
