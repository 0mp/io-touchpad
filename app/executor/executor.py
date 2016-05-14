# -*- coding: utf-8 -*-

"""Executor.

This module executes commands.
"""

import subprocess
import sys
from databox import databox


def execute(command_id):
    """Execute the command related to the command_id.

    This function executes the command it can match with command_id using
    the databox module.

    If the command is found then a child process will be started. It does not
    wait for its termination.

    Args:
        command_id (str): The id which will be matched with a command from
            the databox module.

    Returns:
        True if the related command was found and it tried to launch it. False
        otherwise.
    """
    command_and_arguments = databox.get_command_and_arguments(command_id)
    if command_and_arguments is None:
        print("executor.py: error: the command related to the detected symbol "
              "is missing", file=sys.stderr)
        return False
    else:
        command, arguments = command_and_arguments
        devnull = subprocess.DEVNULL
        subprocess.Popen([command, arguments], stdout=devnull, stderr=devnull)
        return True
