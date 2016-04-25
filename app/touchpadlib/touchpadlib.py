# -*- coding: utf-8 -*-
"""C language touchpadlib library interface for Python."""

import sys
import _thread
import ctypes

LIB_DIRECTORY = "../lib"
TOUCHPADLIB_SHARED_LIBRARY = LIB_DIRECTORY + "/touchpadlib.so"


class TouchpadEvent(ctypes.Structure):
    """The Pythonic version of the struct touchpad_event."""

    _fields_ = [("x", ctypes.c_int32),
                ("y", ctypes.c_int32),
                ("pressure", ctypes.c_int32),
                ("seconds", ctypes.c_long),
                ("useconds", ctypes.c_long)]


class TouchpadSpecification(ctypes.Structure):
    """The Pythonic version of the struct touchpad_specification."""

    _fields_ = [("min_x", ctypes.c_int32),
                ("max_x", ctypes.c_int32),
                ("min_y", ctypes.c_int32),
                ("max_y", ctypes.c_int32),
                ("min_pressure", ctypes.c_int32),
                ("max_pressure", ctypes.c_int32)]


class Touchpadlib:
    """touchpadlib library wrapper.

    Attributes:
        touchpadlib (shared library): The touchpadlib shared library loaded
            with the ctypes module.
        touchpad_event (ctypes.POINTER(TouchpadEvent)): A pointer to the
            class which is a pythonic representation of struct touchpad_event.
        touchpad_specification (ctypes.POINTER(TouchpadSpecification)):
            A pointer to the class which is a pythonic representation of
            struct touchpad_specification.
        touchpad_file_descriptor (ctypes.c_int): A file descriptor to the
            touchpad device.
    """

    touchpadlib = None
    touchpad_event = None
    touchpad_specification = None
    touchpad_file_descriptor = None

    @classmethod
    def get_event(cls):
        """Get the next event from the touchpad.

        Returns:
            Dictionary: It has got the fields of a TouchpadEvent: x, y,
                pressure, seconds and useconds.
        """
        cls.initialize()
        cls.fetch_touchpad_event()
        contents = cls.touchpad_event.contents
        return {'x': contents.x,
                'y': contents.y,
                'pressure': contents.pressure,
                'seconds': contents.seconds,
                'useconds': contents.useconds}

    @classmethod
    def get_specification(cls):
        """Get the specification of a touchpad.

        Returns:
            Dictionary: It has got the fields of a TouchpadSpecification:
                min_x, max_x, min_y, max_y, min_pressure, max_pressure.
        """
        cls.initialize()
        cls.fetch_touchpad_specification()
        contents = cls.touchpad_specification.contents
        return {'min_x': contents.min_x,
                'max_x': contents.max_x,
                'min_y': contents.min_y,
                'max_y': contents.max_y,
                'min_pressure': contents.min_pressure,
                'max_pressure': contents.max_pressure}

    @classmethod
    def initialize(cls):
        """Initialize the static variables."""
        if cls.touchpadlib is None:
            cls.conncect_to_library()
        if cls.touchpad_event is None:
            cls.create_touchpad_event()
        if cls.touchpad_file_descriptor is None:
            cls.initialize_touchpadlib()
        if cls.touchpad_specification is None:
            cls.create_touchpad_specification()

    @classmethod
    def conncect_to_library(cls):
        """Conncet to the touchpadlib."""
        try:
            cls.touchpadlib = \
                    ctypes.cdll.LoadLibrary(TOUCHPADLIB_SHARED_LIBRARY)
        except OSError:
            print("touchpadlib.py: error: no such library as touchpadlib.so",
                  file=sys.stderr)
            cls.interrupt_and_finish()

    @classmethod
    def create_touchpad_event(cls):
        """Get a C language struct touchpad_event from the touchpadlib."""
        cls.touchpadlib.new_event.restype = ctypes.POINTER(TouchpadEvent)
        cls.touchpad_event = cls.touchpadlib.new_event()
        if cls.touchpad_event == 0:
            print("touchpadlib.py: error: cannot allocate memory in "
                  "new_event()", file=sys.stderr)
            cls.interrupt_and_finish()

    @classmethod
    def create_touchpad_specification(cls):
        """Get the touchpad specification struct from the touchpadlib."""
        cls.touchpadlib.new_specification.restype = \
            ctypes.POINTER(TouchpadSpecification)
        cls.touchpad_specification = cls.touchpadlib.new_specification()
        if cls.touchpad_specification == 0:
            print("touchpadlib.py: error: cannot allocate memory in "
                  "new_specification()", file=sys.stderr)
            cls.interrupt_and_finish()

    @classmethod
    def free_touchpad_event(cls):
        """Free the memory allocated for the touchpad_event."""
        cls.touchpadlib.free_event(cls.touchpad_event)

    @classmethod
    def free_touchpad_specification(cls):
        """Free the memory allocated for the touchpad_event."""
        cls.touchpadlib.free_specification(cls.touchpad_specification)

    @classmethod
    def initialize_touchpadlib(cls):
        """Initialize the touchpadlib library.

        Sets the touchpad_file_descriptor.
        """
        cls.touchpadlib.initialize_touchpadlib_usage.restype = ctypes.c_int
        cls.touchpad_file_descriptor = \
            cls.touchpadlib.initialize_touchpadlib_usage()
        if cls.touchpad_file_descriptor == ctypes.c_int(-1):
            print("touchpadlib.py: error: touchpadlib initialize error",
                  file=sys.stderr)
            cls.interrupt_and_finish()

    @classmethod
    def fetch_touchpad_event(cls):
        """Fetch the next event from the touchpad."""
        fetch_event = cls.touchpadlib.fetch_touchpad_event
        file_descriptor = cls.touchpad_file_descriptor
        event = cls.touchpad_event
        if fetch_event(file_descriptor, event) == 1:
            print("touchpadlib.py: error: touchpadlib fetch error",
                  file=sys.stderr)
            cls.interrupt_and_finish()

    @classmethod
    def fetch_touchpad_specification(cls):
        """Fetch the touchpad specification."""
        fetch_specification = cls.touchpadlib.fetch_touchpad_specification
        file_descriptor = cls.touchpad_file_descriptor
        specification = cls.touchpad_specification
        if fetch_specification(file_descriptor, specification) != 0:
            print("touchpadlib.py: error: can't get the touchpad "
                  "specification", file=sys.stderr)
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
            cls.touchpad_event = None
        if cls.touchpad_specification is not None:
            cls.free_touchpad_specification()
            cls.touchpad_specification = None

    @classmethod
    def interrupt_and_finish(cls):
        """Interrupt the main thread and exit.

        Clean the C language data structures.

        Exit thread with the exit code of 1.
        """
        cls.clean()
        _thread.interrupt_main()
        sys.exit(1)
