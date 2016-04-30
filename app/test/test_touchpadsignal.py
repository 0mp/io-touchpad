# -*- coding: utf-8 -*-

"""Tests for touchpad signal """

import pytest

from touchpadsignal import touchpadsignal

# Omit the INIT invokation.
touchpad = touchpadsignal.TouchpadSignal.__new__(touchpadsignal.TouchpadSignal)

def test__combine_seconds_and_useconds():
    seconds = 10
    useconds = 234
    ratio = 0.000001

    result = seconds + useconds * ratio
    method_result = touchpad._combine_seconds_and_useconds(seconds, useconds)

    assert result == method_result

def test_get_x():
    touchpad.x_value = 5
    assert touchpad.get_x() == 5

def test_get_y():
    touchpad.y_value = 6
    assert touchpad.get_y() == 6

def test_get_time():
    touchpad.time = 11
    assert touchpad.get_time() == 11

def test_get_pressure():
    touchpad.pressure = 0.10123
    assert touchpad.get_pressure() == 0.10123

def test_is_proper_signal_of_point():
    for x in range(-1,1):
        for y in range(-1,1):
            touchpad.x_value = x
            touchpad.y_value = y
            ans = True

            if x == -1 or y == -1:
                ans = False

            assert touchpad.is_proper_signal_of_point() == ans



