# -*- coding: utf-8 -*-

"""Tests for classifier"""

import pytest
import operator
from math import fabs

from classifier import classifier as classifier_module

epsilon = 0.00001

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

    assert fabs(classifier.tolerance_distance - 0.5506099238118276) < epsilon

def test_learn():
    """Test learn with existing resource."""
    classifier = classifier_module.Classifier(True,None)
    classifier.training_set_file_path = "learn_dat/training-set.dat"
    classifier.learn(True)
    assert fabs(classifier.tolerance_distance-1271.9887310656133650) < epsilon

def test_classify():
    """Test classify."""
    classifier = classifier_module.Classifier()
    classifier.training_set_file_path = "learn_dat/training-set.dat"
    classifier.training_set = classifier.load_training_set()
    assert classifier.classify(classifier.training_set[0]) == 1

def test__build_paths():
    """Test of the method which build file paths."""
    files1 = ["file1", "file2"]
    userdefined_path = classifier_module.DATA_PATH + classifier_module.USER_DIR

    expected_out_files1 = [operator.add(userdefined_path, file)
                           for file in files1]
    out_files1 = classifier_module.Classifier._build_paths(files=files1,
                                                           system_bitness=None)

    assert len(out_files1) == len(expected_out_files1)

    for file_num in range(len(out_files1)):
        assert out_files1[file_num] == expected_out_files1[file_num]

def test__extend_paths():
    """Test of the tiny helper function which expands a list of strings with a string."""
    file_paths = ["docs/abcd/", "docs/123/"]

    extend_paths = classifier_module.Classifier._extend_paths
    path_element = "u/"
    extended_paths = extend_paths(file_paths, path_element)

    assert len(file_paths) == len(extended_paths)
    for path_num in range(len(file_paths)):
        assert file_paths[path_num] + path_element == extended_paths[path_num]





