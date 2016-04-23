# -*- coding: utf-8 -*-
"""The application thread function.

It receives signals/events from the listener thread using a queue and then
interprets the data.
"""

import time
from signalcollection import signalcollection
from classifier import classifier


def application_thread(queue, learning_mode=False, training_size=0):
    """The application thread function.

    Every iteration of the while loop one signal is read from the queue.
    The signals are sent to the interpreter if there is a longer pause between
    signals or if there is a signal meaning that the hand has been lifted.

    At the moment only the pause is meaningful for us.

    Variables:
        collection (SignalCollection): A collection of signals sent by
            the listener thread. Mostly touchpad events but not necessarily.

    Args:
        queue (Queue): A inter-thread queue to pass signals between the
            listener and the application.
    """
    clsf = classifier.Classifier()
    if learning_mode:
        clsf.reset_training_set(training_size)

    collection = signalcollection.SignalCollection()

    while 1:
        while queue.empty() and collection.is_recent_enough(time.time()):
            pass

        if not collection.is_recent_enough(time.time()):
            send_points_to_interpreter(collection.as_list(), learning_mode, clsf)
            collection.reset()

        if queue.empty():
            continue

        signal = queue.get()

        if signal.is_stop_signal():
            send_points_to_interpreter(collection.as_list(), learning_mode, clsf)
            collection.reset()
        elif signal.is_proper_signal_of_point():
            collection.add_and_maintain(signal)


def send_points_to_interpreter(signal_list, learning_mode, clsf):
    """Interpret the signals from the signal list.

    At the moment the function is not interpreting anything. It just prints the
    first 10 signals from the signal_list. All of them are points.

    Args:
        signal_list (List): List of signals captured from the touchpad.
    """
    if not signal_list:
        return
    print()
    print("New portion of events:")
    counter = 0
    length = len(signal_list)
    for one_signal in signal_list:
        counter += 1
        if counter == 11:
            print("...")
            break
        print("Event #", counter, "\tout of", length, "in the list. x:",
              one_signal.get_x(), "y:", one_signal.get_y())
    print()
    if learning_mode:
        clsf.add_to_training_set(signal_list)
    else:
        item = clsf.classify(signal_list)
        if item is not None:
            # TODO execution
            print("execution")
        else:
            print("not similar")
