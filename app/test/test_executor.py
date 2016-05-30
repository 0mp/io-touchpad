# -*- coding: utf-8 -*-

"""Tests for touchpad signal"""

from executor import executor

def test_execute():
    """Test for execution,

    with unexisting ID and test ID
    """
    default_symbol = '0'
    assert executor.execute('-1') is False
    assert executor.execute(default_symbol) is True
