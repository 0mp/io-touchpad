#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""io-touchpad application.

This application will learn to associate symbols you draw on your touchpad
with certain actions you bound to those symbols.

Your finger's absolute position is fetched using a library called touchpadlib
which is based on evtest (version 1.32).
"""

import argparse
import queue
import threading
import sys

from classifier import classifier as classifier_module
from databox import databox
from terminationhandler import terminationhandler
from threads import application
from threads import listener

MIN_TRAINING_SIZE = 5

ACTIVATE_SUBCOMMAND = 'activate'
ADD_SUBCOMMAND = 'add'
DEACTIVATE_SUBCOMMAND = 'deactivate'
DELETE_SUBCOMMAND = 'delete'
EXPORT_SUBCOMMAND = 'export'
IMPORT_SUBCOMMAND = 'import'
LIST_SUBCOMMAND = 'list'
MODIFY_SUBCOMMAND = 'modify'
REDRAW_SUBCOMMAND = 'redraw'
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

    # The activate subcommand section.
    subparser = subparsers.add_parser(ACTIVATE_SUBCOMMAND,
                                      help='activate selected symbols to be '
                                      'recognisable by the app')
    group = subparser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--all', dest='symbols', action='store_const',
                       const=[], help='activate all the existing symbols')
    group.add_argument('-s', '--select', dest='symbols', nargs="+",
                       help='activate all the selected symbols')

    # The add subcommand section.
    subparser = subparsers.add_parser(ADD_SUBCOMMAND, help='undertake a '
                                      'learning session in order to add '
                                      'a new symbol')
    subparser.add_argument(dest='training_size', metavar='SIZE',
                           help='you will be asked to draw '
                           'a symbol SIZE times; SIZE should be at '
                           'least {0}'.format(MIN_TRAINING_SIZE), type=int)
    subparser.add_argument(dest='symbol_name', metavar='SYMBOL',
                           help='the name of the symbol you want the app to '
                           'learn during the learning session')
    subparser.add_argument(dest='shell_command', metavar='COMMAND',
                           help='the shell command to be triggered when the '
                           'symbol is drawn; for e.g. "touch"')
    subparser.add_argument(dest='shell_command_arguments', metavar='ARGUMENTS',
                           help='the arguments for the shell command to be '
                           'triggered when the symbol is drawn; for '
                           'e.g. "/tmp/touched_file"')

    # The deactivate subcommand section.
    subparser = subparsers.add_parser(DEACTIVATE_SUBCOMMAND, help='deactivate '
                                      'selected symbols; the app will not '
                                      'recognise those symbols; '
                                      'it is be possible to reactivate any of '
                                      'those symbols in the future')
    group = subparser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--all', dest='symbols', action='store_const',
                       const=[], help='deactivate all the existing symbols')
    group.add_argument('-s', '--select', dest='symbols', nargs="+",
                       help='deactivate all the selected symbols')

    # The delete subcommand section.
    subparser = subparsers.add_parser(DELETE_SUBCOMMAND, help='delete '
                                      'selected symbols; the selected symbols '
                                      'will be removed and lost forever')
    group = subparser.add_mutually_exclusive_group(required=True)
    group.add_argument('-a', '--all', dest='symbols', action='store_const',
                       const=[], help='delete all the existing symbols')
    group.add_argument('-s', '--select', dest='symbols', nargs="+",
                       help='delete all the selected symbols')

    # The export subcommand section.
    subparser = subparsers.add_parser(EXPORT_SUBCOMMAND, help='export current '
                                      'settings to a file')
    subparser.add_argument(dest='settings_name',
                           default="settings", metavar='FILE', help='the '
                           'name of the file where the settings will be '
                           'exported')

    # The import subcommand section.
    subparser = subparsers.add_parser(IMPORT_SUBCOMMAND, help='import '
                                      'settings from a file')
    subparser.add_argument(dest='settings_name',
                           default="settings", metavar='FILE', help='the name '
                           'of the file to import settings from')

    # The list subcommand section.
    subparsers.add_parser(LIST_SUBCOMMAND, help='list all the available '
                          'symbols (both activated ' 'and deactivated)')

    # The modify subcommand section.
    subparser = subparsers.add_parser(MODIFY_SUBCOMMAND, help='modify the '
                                      'command assigned to a symbol')
    subparser.add_argument(dest='symbol_name', metavar='SYMBOL', help='the '
                           'name of the symbol the user wants to modify')
    subparser.add_argument(dest='shell_command', metavar='COMMAND', help='the '
                           'new shell command to be triggered when the symbol '
                           'is drawn; for e.g. "touch"')
    subparser.add_argument(dest='shell_command_arguments', metavar='ARGUMENTS',
                           help='the arguments for the shell command to be '
                           'triggered when the symbol is drawn; for '
                           'e.g. "/tmp/touched_file"')

    # The redraw subcommand section.
    subparser = subparsers.add_parser(REDRAW_SUBCOMMAND, help='redraw a '
                                      'symbol while retaining the command '
                                      'related to that symbol; simply, '
                                      'undertake a learning session once '
                                      'again for a specified symbol; the '
                                      'related command, arguments and name '
                                      'will be preserved')
    subparser.add_argument(dest='symbol_name', metavar='SYMBOL', help='the '
                           'name of the symbol you want to draw again')
    subparser.add_argument(dest='training_size', metavar='SIZE',
                           help='you will be asked to redraw '
                           'the selected symbol SIZE times; SIZE should be at '
                           'least {0}'.format(MIN_TRAINING_SIZE), type=int)

    # The repeat subcommand section.
    subparser = subparsers.add_parser(REPEAT_SUBCOMMAND, help='repeat the '
                                      'learning process of all the'
                                      'symbols; use the --symbol option '
                                      'to repeat the process for only one '
                                      'symbol')
    subparser.add_argument('-s', '--symbol', dest='symbol_name',
                           default=None, metavar='SYMBOL', help='repeat '
                           'the learning process for the SYMBOL only; '
                           'if this option is not provided then '
                           'the learning process will be '
                           'repeated for each symbol known to the app')

    # The run subcommand section.
    subparser = subparsers.add_parser(RUN_SUBCOMMAND, help='run the app '
                                      'using either hardcoded or user-defined '
                                      'symbols')
    subparser.add_argument(dest='run_mode', metavar='MODE',
                           choices={RUN_32_MODE, RUN_64_MODE, RUN_USER_MODE},
                           help='set the mode you would like to run the '
                           'app in; use <32> for 32-bit machines, <64> for '
                           '64-bit machines and <user> to use user-defined '
                           'symbols')

    return parser


def _start_threads(learning_mode=None, training_size=None, system_bitness=None,
                   symbol_name=None):
    """Wrap up the threads launching.

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
    condition = threading.Condition()

    # Run both threads.
    listener.start(thread_queue, condition)
    application.application_thread(thread_queue, condition, learning_mode,
                                   training_size, system_bitness, symbol_name)


def _activate(args):
    """Wrap up the activate subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    databox.activate(args.symbols)
    sys.exit(0)


def _add(args):
    """Wrap up the add subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    print('The symbol name is ' + args.symbol_name + '.')
    if args.training_size < MIN_TRAINING_SIZE:
        print('app.py: error: the training size should be at least '
              '{0}'.format(MIN_TRAINING_SIZE), file=sys.stderr)
        sys.exit(1)
    databox.bind_symbol_with_command(args.symbol_name, args.shell_command,
                                     args.shell_command_arguments)
    _start_threads(learning_mode=True, symbol_name=args.symbol_name,
                   training_size=args.training_size)


def _deactivate(args):
    """Wrap up the deactivate subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    databox.deactivate(args.symbols)
    sys.exit(0)


def _delete(args):
    """Wrap up the delete subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    classifier = classifier_module.Classifier(learning_mode=True)
    classifier.delete_symbols(args.symbols)
    databox.delete_symbols(args.symbols)
    sys.exit(0)


def _export_settings(args):
    """Wrap up the export subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    classifier = classifier_module.Classifier(learning_mode=False)
    classifier.export_files(args.settings_name)
    databox.export_settings(args.settings_name)


def _import_settings(args):
    """Wrap up the import subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    classifier = classifier_module.Classifier(learning_mode=True)
    classifier.import_files(args.settings_name)
    databox.import_settings(args.settings_name)


def _list():
    """Wrap up the list subcommand to make main() less complex."""
    databox.print_commands()
    sys.exit(0)


def _modify(args):
    """Wrap up the modify subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    databox.bind_symbol_with_command(args.symbol_name, args.shell_command,
                                     args.shell_command_arguments)
    sys.exit(0)


def _redraw(args):
    """Wrap up the redraw subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    if args.training_size < MIN_TRAINING_SIZE:
        print('app.py: error: the training size should be at least '
              '{0}'.format(MIN_TRAINING_SIZE), file=sys.stderr)
        sys.exit(1)
    _start_threads(learning_mode=True, symbol_name=args.symbol_name,
                   training_size=args.training_size)
    sys.exit(0)


def _repeat(args):
    """Wrap up the repeat subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    print('Repeating the learning process from traning-set file.')
    classifier = classifier_module.Classifier(learning_mode=True)
    classifier.learn(True, args.symbol_name)
    sys.exit(0)


def _run(args):
    """Wrap up the run subcommand to make main() less complex.

    Args:
        args (dict): Parsed command line arguments.
    """
    try:
        system_bitness = int(args.run_mode)
    except ValueError:
        system_bitness = None
    finally:
        _start_threads(system_bitness=system_bitness)


def main():
    """The main function."""
    terminationhandler.setup()

    parser = _get_configured_parser()
    args = parser.parse_args()

    if args.subcommand == ACTIVATE_SUBCOMMAND:
        _activate(args)
    elif args.subcommand == ADD_SUBCOMMAND:
        _add(args)
    elif args.subcommand == DEACTIVATE_SUBCOMMAND:
        _deactivate(args)
    elif args.subcommand == DELETE_SUBCOMMAND:
        _delete(args)
    elif args.subcommand == EXPORT_SUBCOMMAND:
        _export_settings(args)
    elif args.subcommand == IMPORT_SUBCOMMAND:
        _import_settings(args)
    elif args.subcommand == LIST_SUBCOMMAND:
        _list()
    elif args.subcommand == MODIFY_SUBCOMMAND:
        _modify(args)
    elif args.subcommand == REDRAW_SUBCOMMAND:
        _redraw(args)
    elif args.subcommand == REPEAT_SUBCOMMAND:
        _repeat(args)
    elif args.subcommand == RUN_SUBCOMMAND:
        _run(args)

main()
