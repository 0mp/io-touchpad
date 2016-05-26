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

MIN_TRAINING_SIZE = 5

ACTIVATE_SUBCOMMAND = 'activate'
ADD_SUBCOMMAND = 'add'
DEACTIVATE_SUBCOMMAND = 'deactivate'
DELETE_SUBCOMMAND = 'delete'
LIST_SUBCOMMAND = 'list'
MODIFY_SUBCOMMAND = 'modify'
REPEAT_SUBCOMMAND = 'repeat'
RUN_SUBCOMMAND = 'run'
RUN_USER_MODE = 'user'
RUN_32_MODE = '32'
RUN_64_MODE = '64'


def _get_configured_parser():
    """Configure the commandline arguments parser."""
    description = 'Teaching your touchpad magic tricks since 2016.'
    parser = argparse.ArgumentParser(add_help=True, description=description)

    subparsers = parser.add_subparsers(metavar='SUBCOMMAND', dest='subcommand',
                                       help='available subcommands; run '
                                       'SUBCOMMAND -h to show help for '
                                       'the SUBCOMMAND you are interested in')
    subparsers.required = True  # Require SUBCOMMAND to be specified.

    # activate.
    parser_activate = subparsers.add_parser(ACTIVATE_SUBCOMMAND,
                                            help='activate an existing '
                                            'symbol; the app will detect this '
                                            'symbol if the user draws it')
    parser_activate.add_argument(dest='symbol_name', default=None,
                                 metavar='SYMBOL', help='the name of the '
                                 'symbol the user wants to activate')

    # add.
    parser_add = subparsers.add_parser(ADD_SUBCOMMAND, help='undertake a '
                                       'learning session in order to add '
                                       'a new symbol')
    parser_add.add_argument(dest='training_size', default=None,
                            metavar='SIZE', help='you will be asked to draw '
                            'a symbol SIZE times; SIZE should be at '
                            'least ' + str(MIN_TRAINING_SIZE), type=int)
    parser_add.add_argument(dest='symbol_name', default=None,
                            metavar='SYMBOL', help='the name of the symbol '
                            'you want the app to learn during the learning '
                            'session')

    # deactivate.
    parser_deactivate = subparsers.add_parser(DEACTIVATE_SUBCOMMAND,
                                              help='deactivate an existing '
                                              'symbol; the app will not '
                                              'recognise this symbol; it will '
                                              'be possible to activate the '
                                              'symbol back later')
    parser_deactivate.add_argument(dest='symbol_name', default=None,
                                   metavar='SYMBOL', help='the name of the '
                                   'symbol the user wants to deactivate')

    # delete.
    parser_delete = subparsers.add_parser(DELETE_SUBCOMMAND, help='delete an '
                                          'existing symbol')
    parser_delete.add_argument(dest='symbol_name', default=None,
                               metavar='SYMBOL', help='the name of the symbol '
                               'you want to delete from the app')

    # list.
    subparsers.add_parser(LIST_SUBCOMMAND, help='list all the '
                          'available symbols (both activated '
                          'and deactivated)')

    # modify.
    parser_modify = subparsers.add_parser(MODIFY_SUBCOMMAND, help='modify the '
                                          'command assigned to a symbol')
    parser_modify.add_argument(dest='symbol_name', default=None,
                               metavar='SYMBOL', help='the name of the symbol '
                               'the user wants to modify')
    parser_modify.add_argument(dest='new_command', default=None,
                               metavar='COMMAND', help='the new shell command '
                               'to be triggered when the symbol is drawn')

    # repeat.
    parser_repeat = subparsers.add_parser(REPEAT_SUBCOMMAND, help='repeat the '
                                          'learning process of all the'
                                          'symbols; use the --symbol option '
                                          'to repeat the process for only one '
                                          'symbol')
    parser_repeat.add_argument('-s', '--symbol', dest='symbol_name',
                               default=None, metavar='SYMBOL', help='repeat '
                               'the learning process for the SYMBOL only; '
                               'if this option is not provided then '
                               'the learning process will be '
                               'repeated for each symbol known to the app')

    # run.
    parser_run = subparsers.add_parser(RUN_SUBCOMMAND, help='run the app '
                                       'using hardcoded or user-defined '
                                       'symbols')
    parser_run.add_argument(dest='run_mode', metavar='MODE',
                            choices={RUN_32_MODE, RUN_64_MODE, RUN_USER_MODE},
                            help='set the mode you would like to run the '
                            'app in; use <32> for 32-bit machines, <64> for '
                            '64-bit machines and <user> to use user-defined '
                            'symbols')

    return parser


def _start_threads(learning_mode=None, training_size=None, system_bitness=None,
                   symbol_name=None):
    """Wrap the threads launching.

    Args:
        learning_mode (bool): The variable which stores the information if
            the app would like to enter the learning mode or not.
        training_size (int): The number of the learning samples of the symbol
            that the user is asked to draw.
        system_bitness (int): The bitness of the system. The only legal values
            are {None, 32, 64}. If the value is 32 or 64 then set of hardcoded
            symbols (with respect to the provided bitness) will be
            recogniezed instead of the user defined symbols.
        symbol_name (str): The name of the symbol provided by the user with
            a command line option.
    """
    thread_queue = queue.Queue()

    # Run both threads.
    listener.start(thread_queue)
    application.application_thread(thread_queue, learning_mode, training_size,
                                   system_bitness, symbol_name)


def main():
    """The main function."""
    terminationhandler.setup()

    parser = _get_configured_parser()
    args = parser.parse_args()

    if args.subcommand == ACTIVATE_SUBCOMMAND:
        print('app.py: warning: the command line argument "activate SYMBOL" '
              'has not been implemented yet', file=sys.stderr)
        sys.exit(0)
    elif args.subcommand == ADD_SUBCOMMAND:
        print('The symbol name is ' + args.symbol_name + '.')
        if args.training_size < MIN_TRAINING_SIZE:
            print('app.py: error: the training size should be at least %d'
                  % (MIN_TRAINING_SIZE), file=sys.stderr)
            sys.exit(1)
        _start_threads(learning_mode=True, symbol_name=args.symbol_name,
                       training_size=args.training_size)
    elif args.subcommand == DEACTIVATE_SUBCOMMAND:
        print('app.py: warning: the command line argument "deactivate SYMBOL" '
              'has not been implemented yet', file=sys.stderr)
        sys.exit(0)
    elif args.subcommand == DELETE_SUBCOMMAND:
        print('app.py: warning: the command line argument "delete SYMBOL" '
              'has been implemented but for now it only removes symbol'
              'from the classifier', file=sys.stderr)
        classifier = classifier_module.Classifier(True)
        classifier.delete_symbol(args.symbol_name)
        sys.exit(0)
    elif args.subcommand == LIST_SUBCOMMAND:
        print('app.py: warning: the command line argument "list" '
              'has not been implemented yet', file=sys.stderr)
        sys.exit(0)
    elif args.subcommand == MODIFY_SUBCOMMAND:
        print('app.py: warning: the command line argument "modify SYMBOL '
              'COMMAND" has not been implemented yet', file=sys.stderr)
        sys.exit(0)
    elif args.subcommand == REPEAT_SUBCOMMAND:
        print('Repeating the classification within the learning process.')
        classifier = classifier_module.Classifier(True)
        classifier.learn(True, args.symbol_name)
        sys.exit(0)
    elif args.subcommand == RUN_SUBCOMMAND:
        try:
            system_bitness = int(args.run_mode)
        except ValueError:
            system_bitness = None
        finally:
            _start_threads(system_bitness=system_bitness)

main()
