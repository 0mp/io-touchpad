# -*- coding: utf-8 -*-

"""Tests for touchpadlib """

import pytest
import types

from touchpadlib import touchpadlib

def monkey_interrupt_and_finish(self):
    """ Function to moneky patching exiting function.

    It is used to test if something fails"""
    self.ans = -1

#change interupt and finish function for a test
touchpadlib.Touchpadlib.interrupt_and_finish = types.MethodType(monkey_interrupt_and_finish,touchpadlib.Touchpadlib)
touchpadlib.Touchpadlib.ans = 0

tlib = touchpadlib.Touchpadlib

def test_monkey():
    """ Test for changed function

    we check if it changes the value ans"""
    assert tlib.ans == 0
    tlib.interrupt_and_finish()

    assert tlib.ans == -1
    tlib.ans = 0

def test_connect_to_library():
    """ we try to connect to good and bad """
    global TOUCHPADLIB_SHARED_LIBRARY

    #make it wrong
    touchpadlib.TOUCHPADLIB_SHARED_LIBRARY = "sthstupid"
    tlib.connect_to_library()
    assert tlib.ans == -1
    tlib.ans = 0

    #make it right
    touchpadlib.TOUCHPADLIB_SHARED_LIBRARY = "../../lib/touchpadlib.so"
    tlib.connect_to_library()
    assert tlib.ans == 0

    #reset
    tlib.ans = 0
    touchpadlib.TOUCHPADLIB_SHARED_LIBRARY = "../lib/touchpadlib.so"

def test_create_touchpad_event():
    """ test for creating touchpad event """
    tlib.create_touchpad_event()
    assert tlib.ans == 0

def test_initialize_touchpadlib():
    """ test for initialize touchpadlib """
    tlib.initialize_touchpadlib()
    assert tlib.ans == 0


