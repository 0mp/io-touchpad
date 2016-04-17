# -*- coding: utf-8 -*-

"""A collection of signals.

Signal is an item taken out of the queue. It might be something more than a
touchpad event.
"""

MAX_BREAK_BETWEEN_TWO_SIGNALS = 0.3
MAX_NUMBER_OF_SIGNALS_IN_GROUP = 3000
MAX_DURATION_OF_GROUP = 4


class SignalCollection:
    """Collection of signals to interpret.

    Attributes:
        signal_list (list): A list of signals chronologically. Old signals
            are regularly removed.
    """

    def __init__(self):
        """Constructor.

        Simply initializes an empty list.
        """
        self.reset()

    def reset(self):
        """Erase the signal_list."""
        self.signal_list = []

    def add_and_maintain(self, signal):
        """Remove some signals and add a new one.

        Remove signals from the head of the signal_list if they are old
        or if there are too many of them.

        Args:
            signal (TouchpadSignal): A new signal to be added.
        """
        while self.need_to_remove_head(signal):
            self.signal_list.pop(0)
        self.signal_list.append(signal)

    def is_recent_enough(self, current_time):
        """Check if the last signal recent enough.

        Compares the time of the occurance of the tail signal with the
        current time.

        Args:
            current_time (float): The current time which is compared with
                the time of the last signal in the list.
        """
        if self.is_empty():
            return True
        tail_time = self.signal_list[-1].get_time()
        return current_time - tail_time <= MAX_BREAK_BETWEEN_TWO_SIGNALS

    def as_list(self):
        """Return the SignalCollection as a list."""
        return self.signal_list

    def need_to_remove_head(self, signal):
        """Check if there is a need of the list's head removal.

        Return true if the list is too big or if the head is too old.

        Args:
            signal (TouchpadSignal): A signal. Needed to calculate if the head
                is too old or not.
        """
        if self.is_empty():
            return False
        elif self.is_too_big():
            return True
        elif self.is_head_too_old(signal):
            return True
        else:
            return False

    def is_empty(self):
        """True if signal list is empty."""
        return self.signal_list == []

    def is_too_big(self):
        """True if signal list is too big."""
        return len(self.signal_list) >= MAX_NUMBER_OF_SIGNALS_IN_GROUP

    def is_head_too_old(self, signal):
        """True if the head of signal list is too old.

        Args:
            signal (TouchpadSignal): A signal which is compared with the head
                in terms of the time attribute.
        """
        signal_time = signal.get_time()
        head_time = self.signal_list[0].get_time()
        return signal_time - head_time > MAX_DURATION_OF_GROUP
