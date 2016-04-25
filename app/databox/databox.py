# -*- coding: utf-8 -*-
"""Data of the application.

This module stores information about symbols and related commands.

Global variables:
    _BUILTIN_COMMANDS (dict): The id-to-command dictionary of builtin commands.
    _USER_DEFINED_COMMANDS (dict): The id-to-command dictionary of user defined
        commands.
"""


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

_USER_DEFINED_COMMANDS = {}
_BUILTIN_COMMANDS = {
    '0': Command('echo', 'test'),
    '1': Command('x-www-browser', ''),
}


def get_command_and_arguments(command_id):
    """Return the command and arguements related to command_id.

    Args:
        command_id (str): The id of the command.

    Returns:
        The shell command and the arguments related to command_id if the id
        was found. None otherwise.
    """
    if Command.is_builtin(command_id):
        return _BUILTIN_COMMANDS[command_id].get_command_and_argument()
    elif Command.is_user_defined(command_id):
        return _USER_DEFINED_COMMANDS[command_id].get_command_and_argument()
    else:
        return None
