# -*- coding: utf-8 -*-

"""This module is a simple wrapper for signal handling."""

import signal
import sys
from touchpadsignal import touchpadsignal


# pylint: disable=unused-argument
def handler(signum, frame):
    """Free the memory after touchpadlib after SIGINT.

    Args:
        signum (int): Unused.
        frame (int): Unused.
    """
    print("\nClosing the application...")
    touchpadsignal.TouchpadSignal.clean()
    sys.exit(0)


def setup():
    """Set up a SIGINT handler."""
    signal.signal(signal.SIGINT, handler)
