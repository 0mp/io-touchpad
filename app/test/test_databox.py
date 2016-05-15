# -*- coding: utf-8 -*-

"""Tests for touchpad signal"""

import pytest
import os
from databox import databox

command = databox.Command("echo","ok")
DEFAULT_SYMBOL = 'small_a'
DEFAULT_COMMAND = 'x-www-browser'
DEFAULT_ARGUMENTS = ''

# Helper functions

DATABOX_DIR = 'databox'
DATA_DIR = 'databox/data'


# Test section

def test_get_command_and_argument(tmpdir):
    """Test the argument of the command.

    Check if it took good arguments when initialized.
    """
    tmpdir.mkdir(DATABOX_DIR)
    tmpdir.mkdir(DATA_DIR)
    comm, argument = command.get_command_and_argument()

    assert comm == "echo"
    assert argument == "ok"

def test_is_user_defined():
    """Test user defined commands."""
    databox._check_and_load_commands()
    assert command.is_user_defined('not_in') == False

def test_is_builtin():
    """Test the function which tests if the symbol is a builtin.

    False test is invoked with unexisting ID.
    """
    databox._check_and_load_commands()
    assert command.is_builtin(DEFAULT_SYMBOL) == True
    assert command.is_builtin('not_in') == False

def test_get_command_and_arguments(tmpdir):
    """Test the module function for getting commands."""
    tmpdir.mkdir(DATABOX_DIR)
    tmpdir.mkdir(DATA_DIR)
    assert databox.get_command_and_arguments('not_in') == None

    comm, argument = databox.get_command_and_arguments(DEFAULT_SYMBOL)

    assert comm == DEFAULT_COMMAND
    assert argument == DEFAULT_ARGUMENTS

def test__check_and_load_commands(tmpdir):
    tmpdir.mkdir(DATABOX_DIR)
    tmpdir.mkdir(DATA_DIR)
    try:
        databox._check_and_load_commands()
    except:
        assert 0

# def test_bind_symbol_with_command(tmpdir):

