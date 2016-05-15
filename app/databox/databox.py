# -*- coding: utf-8 -*-
"""Data of the application.

This module stores information about symbols and related commands.

Global variables:
    _BUILTIN_COMMANDS (dict): The id-to-command dictionary of builtin commands.
    _USER_DEFINED_COMMANDS (dict): The id-to-command dictionary of user defined
        commands.
"""

import errno
import pickle

class Command(object):
    """A class for storing information about a command.

    This might be an overkill for now but in the future it might allow us
    to expand the functionalities of this class. It will make implementing
    new features to the executor module a lot easier.
    """

    def __init__(self, command, arguments):
        """Constructor.

        Args:
            command (str): The shell command like 'gedit' or 'echo linhchibui'.
            arguments (str): The arguments for the shell command
                like '-la' in 'ls -la'.
        """
        self.command = command
        self.arguments = arguments

    def get_command_and_argument(self):
        """Return the command and the arguments.

        Returns
            The command and the arguments of the object.
        """
        return self.command, self.arguments

    @staticmethod
    def is_user_defined(command_id):
        """Tell if the command is defined in the _USER_DEFINED_COMMANDS.

        Args:
            command_id (str): The id being checked.

        Returns:
            True if is defined in _USER_DEFINED_COMMANDS or False otherwise.
        """
        return command_id in _USER_DEFINED_COMMANDS

    @staticmethod
    def is_builtin(command_id):
        """Tell if the command is defined in the _BUILTIN_COMMANDS.

        Args:
            command_id (str): The id being checked.

        Returns:
            True if is defined in _BUILTIN_COMMANDS or False otherwise.
        """
        return command_id in _BUILTIN_COMMANDS

DATA_PATH = 'databox/data/'
USER_DEFINED_COMMANDS_FILE = 'settings.pickle'

_USER_DEFINED_COMMANDS = None
_BUILTIN_COMMANDS = {
    '0': Command('echo', 'test'),
    'small_a': Command('x-www-browser', ''),
    'large_k': Command('touch', '/tmp/created-by-large-k'),
    'small_gamma': Command('touch', '/tmp/created-by-small_gamma'),
    'small_gamma_with_dot': Command('touch',
                                    '/tmp/created-by-small_gamma_with_dot'),
    'large_sigma': Command('touch', '/tmp/created-by-large_sigma'),
}


def bind_symbol_with_command(symbol, command='touch', command_arguments=None):
    check_and_load_commands()
    if command == 'touch' and command_arguments is None:
        command_arguments = '/tmp/created_by_' + symbol

    _USER_DEFINED_COMMANDS[symbol] = Command(command, command_arguments);
    with open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'wb') as handle:
        pickle.dump(_USER_DEFINED_COMMANDS, handle)

def get_command_and_arguments(command_id):
    """Return the command and arguements related to command_id.

    Args:
        command_id (str): The id of the command.

    Returns:
        The shell command and the arguments related to command_id if the id
        was found. None otherwise.
    """
    check_and_load_commands()
    if Command.is_builtin(command_id):
        return _BUILTIN_COMMANDS[command_id].get_command_and_argument()
    elif Command.is_user_defined(command_id):
        return _USER_DEFINED_COMMANDS[command_id].get_command_and_argument()
    else:
        return None

def check_and_load_commands():
    global _USER_DEFINED_COMMANDS
    if _USER_DEFINED_COMMANDS is None:
        try:
            handle = open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'rb')
        except OSError as e:
            if e.errno == errno.ENOENT:
                _USER_DEFINED_COMMANDS = {}
        else:
                _USER_DEFINED_COMMANDS = pickle.load(handle)

