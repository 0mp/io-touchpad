# -*- coding: utf-8 -*-
"""The Classifier class."""

import math
import operator
import pickle
import sys
import _thread

import numpy as np
from sklearn.neighbors import NearestNeighbors

from classifier import featureextractor

DATA_PATH = 'classifier/data/'

USER_DIR = 'user-defined/'
SYSTEM_BITNESS_32 = 32
SYSTEM_BITNESS_64 = 64
HARDCODED_32BIT_DIR = '32/'
HARDCODED_64BIT_DIR = '64/'

DISTANCE_TOLERANCE_FILE = 'distance-tolerance.dat'
MODEL_FILE = 'nn-model.dat'
TRAINING_SET_FILE = 'training-set.dat'


class Classifier:
    """Class for learning and classifying drawn symbols."""

    def __init__(self, learning_mode=False, system_bitness=None):
        """Constructor. Loads the learning model from files.

        Args:
            learning_mode (bool): Says if we are in the learning mode or not.
            system_bitness (int): The only legal values are {None, 32, 64}.
                If the value is 32 or 64 then set of hardcoded symbols
                (with respect to the provided bitness) will be recogniezed
                instead of the user defined symbols.
        """
        file_names = [DISTANCE_TOLERANCE_FILE, MODEL_FILE, TRAINING_SET_FILE]
        file_paths = Classifier._build_paths(file_names, system_bitness)
        #  (self.distance_tolerance_file_path, self.model_file_path,
        #   self.training_set_file_path) = file_paths

        self.files = {name: path for name, path in zip(file_names, file_paths)}

        if not learning_mode:
            try:
                file_with_model = open(self.files[MODEL_FILE], 'rb')
            except FileNotFoundError:
                print("classifier.py: error: file with the learning model "
                      "doesn't exist; please start the application in the "
                      "learning mode", file=sys.stderr)
                _thread.interrupt_main()
                sys.exit(1)

            self.learning_model = pickle.load(file_with_model)
            file_with_model.close()

            try:
                file_with_tolerance_distance = \
                    open(self.files[DISTANCE_TOLERANCE_FILE], 'r')
            except FileNotFoundError:
                print("classifier.py: error: file with the tolerance distance "
                      "doesn't exist; please start the application in the "
                      "learning mode", file=sys.stderr)
                _thread.interrupt_main()
                sys.exit(1)

            self.tolerance_distance = \
                float(file_with_tolerance_distance.readline())
            file_with_tolerance_distance.close()

        # Variables for learning-mode.
        self.training_size = 0
        self.ultimate_training_size = 0
        self.training_set = []

    def load_training_set(self):
        """Load and return traning symbols from file."""
        try:
            file_with_training = open(self.files[TRAINING_SET_FILE], 'rb')
        except FileNotFoundError:
            print("classifier.py: error: file with training set doesn't "
                  "exist; please start the application in the learning mode",
                  file=sys.stderr)
            _thread.interrupt_main()
            sys.exit(1)

        training_set = pickle.load(file_with_training)
        file_with_training.close()
        return training_set

    def reset_training_set(self, new_training_size):
        """Start the new training set.

        Args:
            new_training_size (int): size of new train-set which have to be
                               given in current learning session.
        """
        self.ultimate_training_size = new_training_size
        self.training_size = 0
        self.training_set = []

    def add_to_training_set(self, signal_list):
        """Add the symbol to training set.

           When all symbols designed for this session are given,
           learning is called.

        Args:
            signal_list (TouchpadSignal list): list of touchpad-signals
            representing the drawn symbol.
        """
        print("training...")
        self.training_set.append(signal_list)
        self.training_size += 1
        print("ok")
        if self.training_size == self.ultimate_training_size:
            self.learn(False)
            _thread.interrupt_main()
            sys.exit(0)
        print()

    def classify(self, signal_list):
        """Classify the symbol to some item.

        Returns item id or None if similirity is too weak.

        Args:
            signal_list (TouchpadSignal list): list of touchpad-signal
            representing the drawn symbol.
        """
        print("classifing...")
        feature_vector = featureextractor.get_features(signal_list)
        distances, _ = self\
            .learning_model.kneighbors(np.array([feature_vector]))
        mean_distance = np.mean(distances[0])
        print(mean_distance)
        if mean_distance < self.tolerance_distance:
            return 1
        else:
            return None

    def compute_tolerance_distance(self, sample):
        """Compute the distance tolerance.

        Returns distance tolerance in the feature vectors space
        below which we find the symbol similar.

        Args:
            sample (list of lists of int): list of feature-vectors,
                                           on which we base on.
        """
        nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree')\
            .fit(sample)
        distances, _ = nbrs.kneighbors(sample)
        print(distances)
        means = []
        for distances_row in distances:
            row = np.delete(distances_row, [0])
            means.append(np.mean(row))
        means.sort()
        critical_index = math.ceil(0.8 * len(means)) - 1
        self.tolerance_distance = means[critical_index] * 1.3
        print("tolerance distance: %.16f" % (self.tolerance_distance))
        file_with_tolerance_distance = \
            open(self.files[DISTANCE_TOLERANCE_FILE], 'w')
        file_with_tolerance_distance.write("%.16f\n"
                                           % (self.tolerance_distance))
        file_with_tolerance_distance.close()

    def learn(self, load_from_file):
        """Learn basing on traing-set.

        Args:
            load_from_file (bool): True - if training has to be load from file,
                          False - new training-set written in self.training_set
                                 that has to be learned and then saved to file.
        """
        print("learning...")
        if not load_from_file:
            file_with_training = open(self.files[TRAINING_SET_FILE], 'wb')
            pickle.dump(self.training_set, file_with_training)
            file_with_training.close()
        training_set = self.load_training_set()
        feature_vectors = []
        for training_element in training_set:
            feature_vectors.append(featureextractor
                                   .get_features(training_element))
        sample = np.array(feature_vectors)
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree')\
            .fit(sample)
        file_with_model = open(self.files[MODEL_FILE], 'wb')
        pickle.dump(nbrs, file_with_model)
        file_with_model.close()
        self.compute_tolerance_distance(sample)

    @staticmethod
    def _build_paths(files, system_bitness):
        """Build paths of the files based on the system bitness.

        Chooses different directories depending on the value of the
        system_bitness. If the bitness is neither 32 nor 64 then the
        USER_DIR directory will be used.

        Args:
            files (list): The names of the files themselves.
            system_bitness (int): The system bitness.
        """
        file_paths = ["" for file in files]

        file_paths = Classifier._extend_paths(file_paths, DATA_PATH)
        if system_bitness == SYSTEM_BITNESS_32:
            file_paths = Classifier._extend_paths(file_paths,
                                                  HARDCODED_32BIT_DIR)
        elif system_bitness == SYSTEM_BITNESS_64:
            file_paths = Classifier._extend_paths(file_paths,
                                                  HARDCODED_64BIT_DIR)
        else:
            file_paths = Classifier._extend_paths(file_paths, USER_DIR)

        file_paths = [operator.add(l, r) for l, r in zip(file_paths, files)]

        return file_paths

    @staticmethod
    def _extend_paths(file_paths, path_element):
        """Extend the file paths with path_element.

        The file paths should end with '/'. Same applies to the path element.

        Args:
            file_paths (list): List of file paths which are going to be
                extended with the path_element.
            path_element (str): The string to be appended to the file_paths.

        Returns:
            file_paths extended with path_element.
        """
        return [operator.add(path, path_element) for path in file_paths]
