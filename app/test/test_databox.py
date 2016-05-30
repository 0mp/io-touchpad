# -*- coding: utf-8 -*-

"""Tests for touchpad signal"""

import pytest
from databox import databox
import shutil

command = databox.Command("echo","ok")
DEFAULT_SYMBOL = 'small_a'
DEFAULT_COMMAND = 'x-www-browser'
DEFAULT_ARGUMENTS = ''

# Helper functions

DATABOX_DIR = 'databox'
DATA_DIR = 'databox/data/'
PRE_TEST_FILE = 'settings.pickle_pretest'

# Setup

def setup_function(test_function):
    """Copy file with prepared data for test"""
    shutil.copy2(DATA_DIR + PRE_TEST_FILE, DATA_DIR + 'settings.pickle')


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
    assert command.is_user_defined('not_in') is False


def test_is_builtin():
    """Test the function which tests if the symbol is a builtin.

    False test is invoked with unexisting ID.
    """
    databox._check_and_load_commands()
    assert command.is_builtin(DEFAULT_SYMBOL) is True
    assert command.is_builtin('not_in') is False


def test_get_command_and_arguments(tmpdir):
    """Test the module function for getting commands."""
    tmpdir.mkdir(DATABOX_DIR)
    tmpdir.mkdir(DATA_DIR)
    assert databox.get_command_and_arguments('not_in') is None

    comm, argument = databox.get_command_and_arguments(DEFAULT_SYMBOL)

    assert comm == DEFAULT_COMMAND
    assert argument == DEFAULT_ARGUMENTS


def test__check_and_load_commands(tmpdir):
    """Test _check_and_load_commands()."""
    tmpdir.mkdir(DATABOX_DIR)
    tmpdir.mkdir(DATA_DIR)
    try:
        databox._check_and_load_commands()
    except:
        assert 0


def test_is_active():
    """Test checking if symbol is active"""
    assert databox.is_active('test')
    assert not databox.is_active('shouldntbe')


def test_single_deactivate():
    """Test deactivating symbol"""
    symbols = ['test']

    databox.deactivate(symbols)
    assert not databox.is_active('test')


def test_single_activate():
    """Test activating symbol"""

    symbols = ['test']

    databox.activate(symbols)
    assert databox.is_active('test')


def test_active_and_deactivate():
    """Test all possibilities with 3 symbols"""

    deactivation_symbols = ['test', 'test2', 'test3']
    for i in range(1, 8):
        databox.deactivate(deactivation_symbols)
        assert not databox.is_active('test')
        assert not databox.is_active('test2')
        assert not databox.is_active('test3')

        to_activate = []
        #activate only selected ones
        if i % 2 == 1:
            to_activate.append('test')

        if i % 4 >= 2:
            to_activate.append('test2')

        if i % 8 >= 4:
            to_activate.append('test3')

        databox.activate(to_activate)
        #independently check if everything is okay
        if i % 2 == 1:
            assert databox.is_active('test')
        else:
            assert not databox.is_active('test')

        if i % 4 >= 2:
            assert  databox.is_active('test2')
        else:
            assert not databox.is_active('test2')

        if i % 8 >= 4:
            assert databox.is_active('test3')
        else:
            assert not databox.is_active('test3')


def test_bind_symbol_with_command():
    """Test binding symbol with command"""

    test = ['test']
    databox.activate(test)
    comm, argument = databox.get_command_and_arguments('test')

    assert comm != 'list'
    assert argument != '--sort'

    databox.bind_symbol_with_command('test', 'list', '--sort')

    comm, argument = databox.get_command_and_arguments('test')

    assert comm == 'list'
    assert argument == '--sort'


def test_delete_symbols():
    """Test deleting symbols"""

    symbols_to_delete = ['test', 'test3']

    databox.delete_symbols(symbols_to_delete)

    temp = ['test']
    databox.activate(temp)
    assert not databox.is_active('test')
    temp = ['test2']
    databox.activate(temp)
    assert databox.is_active('test2')
    temp = ['test3']
    databox.activate(temp)
    assert not databox.is_active('test3')





