# -*- coding: utf-8 -*-
"""A touchpad signal class.

A wraper for the touchpad events from the touchpadlib.

This module communicates with the touchpadlib directly.
"""

import sys
import _thread
from ctypes import cdll, Structure, c_long, POINTER

LIB_DIRECTORY = "../lib"
TOUCHPADLIB_SHARED_LIBRARY = LIB_DIRECTORY + "/touchpadlib.so"


class TouchpadEvent(Structure):
    """The Pythonic version of the touchpad_event struct."""

    _fields_ = [("x", c_long),
                ("y", c_long),
                ("pressure", c_long),
                ("seconds", c_long),
                ("useconds", c_long)]


class TouchpadSignal:
    """Wrapper for the struct touchpad_event from touchpadlib.

    Attributes:
        touchpadlib (a C library): An object which allows to communicate
            with the touchpadlib C library.
        tochpad_event (touchpad_event struct): An struct from the C language
            which stores the fetched event.
        touchpad_file_descriptor: The descriptor to the touchpad device.
    """

    touchpadlib = None
    touchpad_event = None
    touchpad_file_descriptor = None

    def __init__(self):
        """Constructor.

        Firstly, it initialises the needed resources: touchpadlib,
        touchpad_event and touchpad_file_descriptor. Although it will be done
        only once the initialize() function checks whether everythign is ok.

        Secondly, fetch_touchpad_event() captures the latest event from the
        touchpad.
        """
        TouchpadSignal.initialize()
        TouchpadSignal.fetch_touchpad_event()

        touchpadlib = TouchpadSignal.touchpadlib
        touchpad_event = TouchpadSignal.touchpad_event

        self.x_value = touchpadlib.get_x(touchpad_event)
        self.y_value = touchpadlib.get_x(touchpad_event)
        self.pressure = touchpadlib.get_pressure(touchpad_event)
        seconds = touchpadlib.get_seconds(touchpad_event)
        useconds = touchpadlib.get_useconds(touchpad_event)
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

    def is_proper_signal_of_point(self):
        """Check if the signal has both x >= 0 and y >= 0.

        x_value and y_value can be -1 if the event was missing one of them.
        """
        return self.x_value >= 0 and self.y_value >= 0

    @classmethod
    def initialize(cls):
        """Initialize the static variables."""
        if cls.touchpadlib is None:
            cls.conncect_to_library()
        if cls.touchpad_event is None:
            cls.create_touchpad_event()
        if cls.touchpad_file_descriptor is None:
            cls.initialize_touchpadlib()

    @classmethod
    def conncect_to_library(cls):
        """Conncet to the touchpadlib."""
        try:
            cls.touchpadlib = cdll.LoadLibrary(TOUCHPADLIB_SHARED_LIBRARY)
        except OSError:
            print("ERROR: No such library as touchpadlib.so.")
            cls.interrupt_and_finish()

    @classmethod
    def create_touchpad_event(cls):
        """Get a C language struct from the touchpadlib."""
        cls.touchpadlib.new_event.restype = POINTER(TouchpadEvent)
        cls.touchpad_event = cls.touchpadlib.new_event()
        if cls.touchpad_event == 0:
            print("ERROR: Cannot allocate memory in new_event().")
            cls.interrupt_and_finish()

    @classmethod
    def free_touchpad_event(cls):
        """Free the memory allocated for the touchpad_event."""
        cls.touchpadlib.free_event(cls.touchpad_event)

    @classmethod
    def initialize_touchpadlib(cls):
        """Initialize the touchpadlib library.

        Set the touchpad_file_descriptor.
        """
        cls.touchpad_file_descriptor = \
            cls.touchpadlib.initialize_touchpadlib_usage()
        if cls.touchpad_file_descriptor == -1:
            print("ERROR: touchpadlib initialize error.")
            cls.interrupt_and_finish()

    @classmethod
    def fetch_touchpad_event(cls):
        """Fetch the next event from the touchpad."""
        if cls.touchpadlib.fetch_touchpad_event(cls.touchpad_file_descriptor,
                                                cls.touchpad_event) == 1:
            print("ERROR: touchpadlib fetch error.")
            cls.interrupt_and_finish()

    @classmethod
    def clean(cls):
        """Securily free the touchpad_event.

        A function for other modules to be able to free the touchpad_event
        in an emergency.

        terminationhandler uses it.
        """
        if cls.touchpad_event is not None:
            cls.free_touchpad_event()

    @classmethod
    def interrupt_and_finish(cls):
        """Interrupt the main thread and exit.

        Clean the C language data structures.

        Exit thread with the exit code of 1.
        """
        cls.clean()
        _thread.interrupt_main()
        sys.exit(1)

    @staticmethod
    def combine_seconds_and_useconds(seconds, useconds):
        """Combine seconds and miliseconds into one variable."""
        return seconds + 0.000001 * useconds
