# -*- coding: utf-8 -*-

"""Tests for signal collection"""

import pytest

from classifier import featureextractor
from math import pi

# Basic points
POINTA = [2.0, 4.0]
POINTB = [1.0, 2.0]

def test_center_of_line():
    """Simple center of line test
    
    """
    assert featureextractor.center_of_line(POINTA, POINTB) == (1.5, 3.0)

def test_length_of_line():
    """Simple length of line test

    """
    assert featureextractor.length_of_line(POINTA, POINTB) == 2.23606797749979

# Class for mocking signal
class Signal_test:
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def get_x(self): 
        return self.x
    
    def get_y(self):
        return self.y

    def is_proper_signal_of_point(self):
        return True

    def is_raising_finger_signal(self):
        return False

SIGNALA = Signal_test(1.0, 1.0)
SIGNALB = Signal_test(2.0, 2.0)
SIGNAL_LIST_TEST = [SIGNALA, SIGNALB]

def test_calculate_border_points():
    """Test the border point of signal test list

    """
    assert featureextractor.calculate_border_points(SIGNAL_LIST_TEST) == (1.0,1.0,2.0,2.0)

def test_calculate_center_of_mass_and_length():
    """Test for calculate center of mass and length

    """
    assert featureextractor.calculate_center_of_mass_and_length(SIGNAL_LIST_TEST) == ((1.5, 1.5), 1.4142135623730951)

def test_ratio_point_of_line():
    """Test for ratio point of line function
    
    """
    assert featureextractor.ratio_point_of_line(POINTA, POINTB, 0.5) == (1.5, 3.0)

def test_scale_point():
    """Test for scale point function

    """
    assert featureextractor.scale_point(POINTA, 0.0, 0.0, 0.0, 0.0, POINTB) == (0, 0)
    assert featureextractor.scale_point(POINTA, 0.0, 0.0, 10.0, 10.0, POINTB) == (100.0, 200.0)

def test_get_angle_between_line_and_xaxis():
    """Test for angle between line and xaxis

    """
    PA = [2, 2]
    PB = [-6, -6]
    assert featureextractor.get_angle_between_line_and_xaxis(PA,PB) == pi/4


def test_join_features():
    """Test for joining features

    """
    P1 = [3, 4]
    P2 = [5, 6]
    P3 = [7, 8]
    PLIST = [P1, P2, P3]
    F1 = [3.3, 2.2]
    C1 = [5.6, 45, 12]

    assert len(featureextractor.join_features(PLIST, F1, C1)) == 11

