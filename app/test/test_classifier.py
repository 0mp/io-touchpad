# -*- coding: utf-8 -*-

"""Tests for classifier"""

import pytest

from classifier import classifier as classifier_module



def test_reset_training_set():
    """Test for load training set."""
    classifier = classifier_module.Classifier(True,None)
    classifier.reset_training_set(117)
    assert classifier.training_set == []
    assert classifier.training_size == 0
    assert classifier.ultimate_training_size == 117


def test_compute_tolerance_distance():
    """Test for computing distance.

    We put some list of list of features to calculate 
    fixed distance, and check if it's same.
    """
    classifier = classifier_module.Classifier(True,None)
    L1 = [11.2, 41.43, 1.33]
    L2 = [10.9, 41.45, 1.34]
    L3 = [12.0, 41.4412, 1.001]
    L4 = [11.3, 41.15, 1.12]
    L5 = [11.223, 41.0, 1.31]
    AL = [L1, L2, L3, L4, L5]
    classifier.compute_tolerance_distance(AL)

    assert classifier.tolerance_distance == 0.5506099238118276

def test_learn():
    """Test learn with existing resource."""
    classifier = classifier_module.Classifier(True,None)
    classifier.training_set_file_path = "learn_dat/training-set.dat"
    classifier.learn(True)
    assert classifier.tolerance_distance == 1271.9887310656133650

def test_classify():
    """Test classify."""
    classifier = classifier_module.Classifier()
    classifier.training_set_file_path = "learn_dat/training-set.dat"
    classifier.training_set = classifier.load_training_set()
    assert classifier.classify(classifier.training_set[0]) == 1
