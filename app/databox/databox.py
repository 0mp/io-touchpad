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

    def __init__(self, command, arguments, status='active'):
        """Constructor.
        Args:
            command (str): The shell command like 'gedit' or 'echo linhchibui'.
            arguments (str): The arguments for the shell command
                like '-la' in 'ls -la'.
        """
        self.command = command
        self.arguments = arguments
        self.status = status

    def to_str(self):
        """Return string containing information about the command."""
        return self.command + ' ' + self.arguments + ' ' + self.status

    def is_active(self):
        """Return True if the command is active. False othervise"""
        return self.status == 'active'

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


def is_active(symbol):
    """Check if given symbol is active.

    if there is no such symbol in database, returns false"""
    global _USER_DEFINED_COMMANDS
    _check_and_load_commands()

    if Command.is_builtin(symbol):
        command = _BUILTIN_COMMANDS[symbol]
    elif Command.is_user_defined(symbol):
        command = _USER_DEFINED_COMMANDS[symbol]
    else:
        return None

    return command.is_active()


def print_commands():
    global _USER_DEFINED_COMMANDS
    _check_and_load_commands()
    for sym in _USER_DEFINED_COMMANDS:
        command = _USER_DEFINED_COMMANDS[sym]
        print(sym, command.to_str())


def _set_status(symbols, status):
    global _USER_DEFINED_COMMANDS
    _check_and_load_commands()
    if not symbols:
        symbols = _USER_DEFINED_COMMANDS
    for symbol in symbols:
        if symbol in _USER_DEFINED_COMMANDS:
            _USER_DEFINED_COMMANDS[symbol].status = status
        else:
            print('warning: symbol', symbol, 'is not present in databox')
    with open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'wb') as handle:
        pickle.dump(_USER_DEFINED_COMMANDS, handle)


def activate(symbols):
    """Change the status of given symbols on active."""
    _set_status(symbols, 'active')


def deactivate(symbols):
    """Change the status of given symbols on inactive."""
    _set_status(symbols, 'inactive')


def delete_symbols(symbols):
    """Delete commands related to given symbols.
    if symbols is empty list, then remove all symbols.
    Args:
        symbols (list of str): Symbols to remove names.
    """
    global _USER_DEFINED_COMMANDS
    if not symbols:
        print('removing all symbols from databox...')
        _USER_DEFINED_COMMANDS = {}
    elif symbols: 
        _check_and_load_commands()
        for symbol in symbols:
            print('removing symbol', symbol, 'from databox')
            if symbol in _USER_DEFINED_COMMANDS:
                del _USER_DEFINED_COMMANDS[symbol]
            else:
                print('warning: symbol', symbol, 'is not present in databox')
    with open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'wb') as handle:
        pickle.dump(_USER_DEFINED_COMMANDS, handle)
         

def bind_symbol_with_command(symbol, command='touch', command_arguments=None):
    """Bind the symbol's name with the provided command.
    The default command is touch and the command_arguments will be set to
    '/tmp/created_by_' + symbol.
    Args:
        symbol (str): The symbol's name.
        command (str): The shell command to be bound with the symbol.
        command_arguments (str): The arguments for the command.
    """
    global _USER_DEFINED_COMMANDS
    _check_and_load_commands()
    if command == 'touch' and command_arguments is None:
        command_arguments = '/tmp/created_by_' + symbol

    _USER_DEFINED_COMMANDS[symbol] = Command(command, command_arguments)
    with open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'wb') as handle:
        pickle.dump(_USER_DEFINED_COMMANDS, handle)


def get_command_and_arguments(command_id):
    """Return the command and arguements related to command_id.
    Args:
        command_id (str): The id of the command.
    Returns:
        The shell command and the arguments related to command_id if the id
        was found and is active. None otherwise.
    """
    _check_and_load_commands()
    if Command.is_builtin(command_id):
        command = _BUILTIN_COMMANDS[command_id]
    elif Command.is_user_defined(command_id):
        command = _USER_DEFINED_COMMANDS[command_id]
    else:
        command = None
    if command is not None and command.is_active():
        return command.get_command_and_argument()
    else:
        return None


def _check_and_load_commands():
    """Check whether the user-defined commands've been loaded and load them.
    It uses the global statement but it is unavoidable in the current
    architecture of the databox module.
    """
    global _USER_DEFINED_COMMANDS
    if _USER_DEFINED_COMMANDS is None:
        try:
            handle = open(DATA_PATH + USER_DEFINED_COMMANDS_FILE, 'rb')
        except OSError as file_not_found_error:
            if file_not_found_error.errno == errno.ENOENT:
                _USER_DEFINED_COMMANDS = {}
        else:
            _USER_DEFINED_COMMANDS = pickle.load(handle)
