# -*- coding: utf-8 -*-

"""Tests for touchpad signal"""

from executor import executor

def test_execute():
    """Test for execution,

    with unexisting ID and test ID
    """
    default_symbol = 'small_a'
    assert executor.execute('-1') == False
    assert executor.execute('small_a') == True
