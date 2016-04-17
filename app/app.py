#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""io-touchpad application.

This application will learn to associate symbols you draw on your touchpad
with certain actions you bound to those symbols.

Your finger's absolute position is fetched using a library called touchpadlib
which is based on evtest (version 1.32).

Source code encoding: UTF-8
"""

import queue

from terminationhandler import terminationhandler
from threads import application
from threads import listener


def main():
    """The main function."""
    # SIGINT signal handler.
    terminationhandler.setup()

    thread_queue = queue.Queue()

    # Run both threads.
    print("\nUse your touchpad as usual. Have a nice day!")

    listener.start(thread_queue)
    application.application_thread(thread_queue)

main()
