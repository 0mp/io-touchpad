# -*- coding: utf-8 -*-
"""Classifier class."""

import _thread
import pickle
from sklearn.neighbors import NearestNeighbors
import numpy as np
import sys
import math
from classifier import normalizer


class Classifier:

    """Class for learning and classifying drawn symbols."""

    def __init__(self):
        """Constructor. Loads learning model from file."""
        file_with_model = open('nn-model', 'rb')
        self.learning_model = pickle.load(file_with_model)
        file_with_model.close()
        file_with_tolerance_distance = open('tolerance_distance', 'r')
        self.tolerance_distance = float(file_with_tolerance_distance.readline())
        file_with_tolerance_distance.close()
        self.training_size = 0
        self.ultimate_training_size = 0
        self.file_with_sizes = None
        self.file_with_signals = None

    def load_training_set(self):
        """Load traning symbols from file."""
        self.file_with_sizes = open('drawings-sizes', 'r')
        self.file_with_signals = open('drawings-signals', 'rb')
        number_of_symbols = int(self.file_with_sizes.readline())
        training_set = []
        for _ in range(0, number_of_symbols):
            symbol_size = int(self.file_with_sizes.readline())
            signals = []
            for _ in range(0, symbol_size):
                touchpad_signal = pickle.load(self.file_with_signals)
                signals.append(touchpad_signal)
            training_set.append(signals)
        self.file_with_sizes.close()
        self.file_with_signals.close()
        return training_set

    def reset_training_set(self, new_training_size):
        """Start the new training set."""
        self.ultimate_training_size = new_training_size
        self.training_size = 0
        self.file_with_sizes = open('drawings-sizes', 'w')
        self.file_with_sizes.write("%d\n" % (new_training_size))
        self.file_with_signals = open('drawings-signals', 'wb')

    def add_to_training_set(self, signal_list):
        """Add the symbol to training set."""
        print("training...")
        self.file_with_sizes.write("%d\n" % (len(signal_list)))
        for element in signal_list:
            pickle.dump(element, self.file_with_signals)
        self.training_size += 1
        print("ok")
        if self.training_size == self.ultimate_training_size:
            self.file_with_sizes.close()
            self.file_with_signals.close()
            self.learn()
            _thread.interrupt_main()
            sys.exit(0)
        print()

    def calculate_feature_vector(self, signal_list):
        """Calculate vector of features for given symbol."""
        # temporal stupid features:
        length = len(signal_list)
        feature_vector = []
        for i in range(0, 30):
            index = int((length * i) / 30)
            feature_vector.append(signal_list[index].get_x())
            feature_vector.append(signal_list[index].get_y())
        return feature_vector

    def classify(self, signal_list):
        """Classify the symbol to some item id or return None if similirity is to weak."""
        print("classyfing...")
        feature_vector = normalizer.get_features(signal_list)
        # TODO normalizing features by variance or spread
        distances, _ = self.learning_model.kneighbors(np.array([feature_vector]))
        mean_distance = np.mean(distances[0])
        if mean_distance < self.tolerance_distance:
            return 1
        else:
            return None

    def compute_tolerance_distance(self, sample):
        """Compute the distance in the feature vectors space below which we find the symbol similar."""
        nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree').fit(sample)
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
        file_with_tolerance_distance = open('tolerance_distance', 'w')
        file_with_tolerance_distance.write("%.16f\n" % (self.tolerance_distance))
        file_with_tolerance_distance.close()

    def learn(self):
        """Load training symbols and learn."""
        print("learning...")
        training_set = self.load_training_set()
        feature_vectors = []
        for training_element in training_set:
            # TODO normalizing fetaures by variance or spread
            feature_vectors.append(normalizer.get_features(training_element))
        sample = np.array(feature_vectors)
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(sample)
        file_with_model = open('nn-model', 'wb')
        pickle.dump(nbrs, file_with_model)
        file_with_model.close()
        self.compute_tolerance_distance(sample)
