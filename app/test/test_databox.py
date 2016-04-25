# -*- coding: utf-8 -*-

"""Tests for touchpad signal"""

import pytest

from databox import databox

command = databox.Command("echo","ok")

# Test section

def test_get_command_and_argument():
    """Test for argument function, if it took good arguments when initialized

    """
    comm, argument = command.get_command_and_argument()
    
    assert comm == "echo"
    assert argument == "ok"

def test_is_user_defined():
    """Test for user defined functions

    """
    assert command.is_user_defined('not_in') == False

def test_is_builtin():
    """Test for built in function,

    False test is invoked with unexisting ID.
    """
    assert command.is_builtin('1') == True
    assert command.is_builtin('not_in') == False

def test_get_command_and_arguments():
    """Test for module function for getting commands

    """
    assert databox.get_command_and_arguments('not_in') == None

    comm, argument = databox.get_command_and_arguments('1')
    
    assert comm == 'x-www-browser'
    assert argument == ''
