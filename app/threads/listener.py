# -*- coding: utf-8 -*-
"""The listener thread.

This module has code of the listener thread.

The job of this thread is to fetch events from the touchpad using touchpadlib
and put them onto a queue so that the other thread could interpret the events.
"""

import _thread

from touchpadsignal import touchpadsignal


def start(queue):
    """Start the listener thread.

    Args:
        queue, (Queue): An inter-thread communication queue for threads.
    """
    _thread.start_new_thread(listener_thread, (queue,))


def listener_thread(queue):
    """The main function of the listener thread.

    The listener thread fetches events from the touchpad and pushes them onto
    the queue.

    Args:
        queue, (Queue): An inter-thread communication queue for threads.
    """
    while 1:
        touchpad_signal = touchpadsignal.TouchpadSignal()
        queue.put(touchpad_signal, True)
