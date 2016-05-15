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
from classifier import classifier as classifier_module

MIN_TRAINING_SIZE = 3

RUN_SUBCOMMAND = 'run'
RUN_USER_MODE = 'user'
RUN_32_MODE = '32'
RUN_64_MODE = '64'
LEARN_SUBCOMMAND = 'learn'
REPEAT_SUBCOMMAND = 'repeat'

def _get_configured_parser():
    """Configure the commandline arguments parser."""
    description = 'Teaching your touchpad magic tricks since 2016.'
    parser = argparse.ArgumentParser(add_help=True, description=description)

    subparsers = parser.add_subparsers(metavar='SUBCOMMAND', dest='subcommand',
            help='available subcommands; run <subcommand> -h to show help for '
            'the <subcommand>')
    subparsers.required = True # Require SUBCOMMAND to be specified.

    parser_run = subparsers.add_parser(RUN_SUBCOMMAND, help='a subcommand to '
            'run the app using hardcoded or user-defined symbols')
    parser_run.add_argument(dest='run_mode', metavar='MODE',
            choices={RUN_32_MODE, RUN_64_MODE, RUN_USER_MODE},
            help='set the mode you would like to run the app in; use <32> '
            'for 32-bit machines, <64> for 64-bit machines and <user> to use '
            'user-defined symbols')

    parser_learn = subparsers.add_parser(LEARN_SUBCOMMAND,
            help='a subcommand to undertake a learning session')
    parser_learn.add_argument(dest='training_size', default=None,
            metavar='SIZE', help='you will be asked to draw a symbol SIZE '
            'times; SIZE should be at least ' + str(MIN_TRAINING_SIZE),
            type=int)
    parser_learn.add_argument(dest='symbol_name', default=None,
            metavar='SYMBOL', help='the name of the symbol ' 'you want the '
            'app to learn during the learning session')

    parser_repeat = subparsers.add_parser(REPEAT_SUBCOMMAND,
            help='a subcommand to repeat the classification process of all '
            'symbols; use the --symbol option to repeat the process for '
            'only one symbol.')
    parser_repeat.add_argument('-s', '--symbol', dest='symbol_name',
            default=None, metavar='SYMBOL', help='the name of the symbol on '
            'which you want to repeat the classification process')

    return parser


def main():
    """The main function."""
    # SIGINT signal handler.
    terminationhandler.setup()

    parser = _get_configured_parser()
    args = parser.parse_args()

    training_size = None
    system_bitness = None
    learning_mode = None
    symbol_name = None

    #  if args.repeat_classification:
    if args.subcommand == REPEAT_SUBCOMMAND:
        print('Repeating the classification within the learning process.')
        classifier = classifier_module.Classifier(True)
        classifier.learn(True, args.symbol_name)
        sys.exit(0)

    #  if learning_mode:
    if args.subcommand == LEARN_SUBCOMMAND:
        print('The symbol name is ' + args.symbol_name + '.')
        if training_size < 3:
            print('app.py: error: the training size should be at least %d'
                  % (MIN_TRAINING_SIZE), file=sys.stderr)
            sys.exit(1)
        training = args.training_size
        symbol_name = args.symbol_name
        learning_mode = True

    if args.subcommand == RUN_SUBCOMMAND:
        if args.run_mode != RUN_USER_MODE:
            system_bitness = int(args.run_mode)


    thread_queue = queue.Queue()

    # Run both threads.
    listener.start(thread_queue)
    application.application_thread(thread_queue, learning_mode, training_size,
                                   system_bitness, symbol_name)

main()
