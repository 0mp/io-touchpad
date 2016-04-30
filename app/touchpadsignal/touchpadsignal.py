# -*- coding: utf-8 -*-
"""A touchpad signal class.

A wraper for the touchpad events from the touchpadlib.

This module communicates with the touchpadlib directly.
"""

from touchpadlib import touchpadlib


class TouchpadSignal:
    """Wrapper for the struct touchpad_event from touchpadlib.

    Attributes:
        touchpadlib (a C library): An object which allows to communicate
            with the touchpadlib C library.
        tochpad_event (touchpad_event struct): An struct from the C language
            which stores the fetched event.
        touchpad_file_descriptor: The descriptor to the touchpad device.
    """

    def __init__(self):
        """Constructor.

        Firstly, it initialises the needed resources: touchpadlib,
        touchpad_event and touchpad_file_descriptor. Although it will be done
        only once the initialize() function checks whether everythign is ok.

        Secondly, fetch_touchpad_event() captures the latest event from the
        touchpad.
        """
        #  TouchpadSignal.initialize()
        #  TouchpadSignal.fetch_touchpad_event()

        #  touchpadlib = TouchpadSignal.touchpadlib
        #  touchpad_event = TouchpadSignal.touchpad_event

        touchpad_event = touchpadlib.Touchpadlib.get_event()

        self.x_value = touchpad_event['x']
        self.y_value = touchpad_event['y']
        self.pressure = touchpad_event['pressure']
        seconds = touchpad_event['seconds']
        useconds = touchpad_event['useconds']
        self.time = self.combine_seconds_and_useconds(seconds, useconds)

    def get_x(self):
        """Get the x value."""
        return self.x_value

    def get_y(self):
        """Get the y value."""
        return self.y_value

    def get_time(self):
        """Get the time."""
        return self.time

    def get_pressure(self):
        """Get the pressure."""
        return self.pressure

    def is_stop_signal(self):
        """Check if the signal is a stop signal."""
        #  @TODO Implement this method. Not in the first iteration.
        #  We could set condition pressure<1 but for now the only reason
        #  to end the group is a long break between signals.
        return False

    def is_raising_finger_signal(self):
        """ Check if the signal has pressure equal to 0 """
        return self.pressure == 0

    def is_proper_signal_of_point(self):
        """Check if the signal has both x >= 0 and y >= 0.

        x_value and y_value can be -1 if the event was missing one of them.
        """
        return self.x_value >= 0 and self.y_value >= 0

    @staticmethod
    def combine_seconds_and_useconds(seconds, useconds):
        """Combine seconds and miliseconds into one variable."""
        return seconds + 0.000001 * useconds
