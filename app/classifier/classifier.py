# -*- coding: utf-8 -*-
"""The Classifier class."""

import math
import operator
import pickle
import sys
import _thread
import os

import numpy as np
from sklearn.neighbors import NearestNeighbors

from sklearn.neighbors import KNeighborsClassifier

from classifier import featureextractor

from string import Template

DATA_PATH = 'classifier/data/'
TEST_LOCATION = 'classifier/data/'

USER_DIR = 'user-defined/'
SYSTEM_BITNESS_32 = 32
SYSTEM_BITNESS_64 = 64
HARDCODED_32BIT_DIR = '32/'
HARDCODED_64BIT_DIR = '64/'

EXPORT_DIR = 'exports/'

DISTANCE_TOLERANCE_FILE = 'distance-tolerance$sym.dat'
EXPORT_SAVING_FILE = 'exports/$sym'
MODEL_FILE = 'nn-model$sym.dat'
TRAINING_SET_FILE = 'training-set$sym.dat'
SYMBOL_LIST_FILE = 'symbol-list.dat'
INACTIVE_SYMBOLS_FILE = 'inactive-symbols.dat'

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
        file_names = [DISTANCE_TOLERANCE_FILE, MODEL_FILE,
                      TRAINING_SET_FILE, SYMBOL_LIST_FILE,
                      EXPORT_SAVING_FILE, INACTIVE_SYMBOLS_FILE]

        file_paths = Classifier._build_paths(file_names, system_bitness)
        #  (self.distance_tolerance_file_path, self.model_file_path,
        #   self.training_set_file_path) = file_paths
        self.files = {name: path for name, path in zip(file_names, file_paths)}

        # Symbol list loading.
        try:

            with open(self.files[SYMBOL_LIST_FILE], 'rb') as handle:
                self.symbol_list = pickle.load(handle)

        except FileNotFoundError:
            self.symbol_list = []

        # Loading classifying models and distance tolerances
        self.learning_models = []
        self.tolerance_distances = []
        self.symbol_list.append("")
        for symbol in self.symbol_list:
            try:

                model_path = Classifier.\
                    _get_file_path(self.files[MODEL_FILE], symbol)

                with open(model_path, 'rb') as handle:
                    self.learning_models.append(pickle.load(handle))

            except FileNotFoundError:
                if learning_mode:
                    self.learning_models.append({})
                else:
                    print("classifier.py: error: file with the learning model "
                          "doesn't exist; please start the application in the "
                          "learning mode", file=sys.stderr)
                    print(model_path)
                    _thread.interrupt_main()
                    sys.exit(1)

            if symbol != "":
                try:

                    tolerance_distance_path = \
                        Classifier._get_file_path(
                            self.files[DISTANCE_TOLERANCE_FILE], symbol)

                    with open(tolerance_distance_path, 'r') as handle:
                        self.tolerance_distances.\
                            append(float(handle.readline()))

                except FileNotFoundError:
                    if learning_mode:
                        self.tolerance_distances.append(0)
                    else:
                        print("classifier.py: error: file with the tolerance "
                              "distance doesn't exist; please start the "
                              "application in the learning mode",
                              file=sys.stderr)
                        print(tolerance_distance_path)
                        _thread.interrupt_main()
                        sys.exit(1)

        self.learning_models = \
            {sym: mod for sym, mod in
             zip(self.symbol_list, self.learning_models)}
        self.tolerance_distances = \
            {sym: dist for sym, dist in
             zip(self.symbol_list, self.tolerance_distances)}
        self.symbol_list.pop()

        # Variables for learning-mode.
        self.training_size = 0
        self.ultimate_training_size = 0
        self.training_set = []
        self.symbol_name = None

    def _load_training_set(self, symbol, file_not_found_ignore=False):
        """Load and return traning symbols from file."""
        try:
            training_path = Classifier.\
                _get_file_path(self.files[TRAINING_SET_FILE], symbol)

            with open(training_path, 'rb') as handle:
                training_set = pickle.load(handle)

        except FileNotFoundError:
            if file_not_found_ignore:
                return None
            print("classifier.py: error: file with training set doesn't "
                  "exist; please start the application in the learning mode",
                  file=sys.stderr)
            print(training_path)
            _thread.interrupt_main()
            sys.exit(1)

        return training_set

    def export_files(self, settings_name):
        """Export saved settings to file.

        Args:
            settings_name (str): The id of the saved settings.
        """
        print('exporting in classifier')
        box = []
        box.append(self.symbol_list)
        inactive_symbols = self._get_inactive_symbols()
        box.append(inactive_symbols)
        box.append(self.learning_models[""])
        for symbol in self.symbol_list:
            box.append(self.learning_models[symbol])
            box.append(self.tolerance_distances[symbol])
            training_set = self._load_training_set(symbol, True)
            box.append(training_set)
        export_path = Classifier.\
            _get_file_path(self.files[EXPORT_SAVING_FILE], settings_name)
        if not os.path.exists(DATA_PATH + USER_DIR + EXPORT_DIR):
            os.makedirs(DATA_PATH + USER_DIR + EXPORT_DIR)
        file_with_export = open(export_path, 'wb')
        pickle.dump(box, file_with_export)
        file_with_export.close()

    def import_files(self, settings_name):
        """Import saved settings from file.

        Args:
            settings_name (str): The id of the saved settings.
        """
        print('importing in classifier')
        try:
            export_path = Classifier.\
                _get_file_path(self.files[EXPORT_SAVING_FILE], settings_name)
            file_with_export = open(export_path, 'rb')
            box = pickle.load(file_with_export)
            file_with_export.close()
        except FileNotFoundError:
            print("name of settings not found in classifier database")
            _thread.interrupt_main()
            sys.exit(1)
        symbol_list = box.pop(0)
        file_with_symbols = \
            open(self.files[SYMBOL_LIST_FILE], 'wb')
        pickle.dump(symbol_list, file_with_symbols)
        file_with_symbols.close()

        inactive_symbols = box.pop(0)
        file_with_inactive_symbols = \
            open(self.files[INACTIVE_SYMBOLS_FILE], 'wb')
        pickle.dump(inactive_symbols, file_with_inactive_symbols)
        file_with_inactive_symbols.close()

        general_model = box.pop(0)
        file_with_model = \
            open(Classifier._get_file_path(self.files[MODEL_FILE], ""),
                 'wb')
        pickle.dump(general_model, file_with_model)
        file_with_model.close()

        for symbol in symbol_list:
            learning_model = box.pop(0)
            file_with_model = \
                open(Classifier._get_file_path(self.files[MODEL_FILE], symbol),
                     'wb')
            pickle.dump(learning_model, file_with_model)
            file_with_model.close()

            tolerance_distance = box.pop(0)
            tolerance_distance_path = Classifier._get_file_path(
                self.files[DISTANCE_TOLERANCE_FILE], symbol)
            file_with_tolerance_distance = \
                open(tolerance_distance_path, 'w')
            file_with_tolerance_distance.write("%.16f\n"
                                               % (tolerance_distance))
            file_with_tolerance_distance.close()

            training_set = box.pop(0)
            file_with_training_path = \
                Classifier._get_file_path(self.files[TRAINING_SET_FILE], symbol)
            file_with_training = \
                open(file_with_training_path, 'wb')
            pickle.dump(training_set, file_with_training)
            file_with_training.close()

    def reset_training_set(self, new_training_size, symbol_name):
        """Start the new training set.

        Args:
            new_training_size (int): size of new train-set which have to be
                               given in current learning session.
            symbol_name    (String): name of the symbol which training set
                               is being resetted.
        """
        self.ultimate_training_size = new_training_size
        self.training_size = 0
        self.training_set = []
        self.symbol_name = symbol_name

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
            self.learn(False, self.symbol_name)
            _thread.interrupt_main()
            sys.exit(0)
        print()

    def classify(self, signal_list):
        """Classify the symbol to some an item.

        Args:
            signal_list (TouchpadSignal list): the list of the signals fetched
                from a touchpad representing the drawn symbol.

        Returns:
            The name of the symbol (such as "small_a" for a or "large_k for K
            if similarity has been found. None otherwise.
        """
        print("classifying...")
        feature_vector = featureextractor.get_features(signal_list)
        models = self.learning_models
        if models[""] == {}:
            return None
        symbol_candidate = models[""].predict([feature_vector])[0]
        if models[symbol_candidate] == {}:
            return None
        distances, _ = models[symbol_candidate] \
            .kneighbors(np.array([feature_vector]))
        mean_distance = np.mean(distances[0])
        print(mean_distance)
        if mean_distance < self.tolerance_distances[symbol_candidate]:
            print(symbol_candidate)
            return symbol_candidate
        else:
            return None

    def _compute_tolerance_distance(self, sample, symbol):
        """Compute the distance tolerance.

        Computes distance tolerance in the feature vectors space
        below which we find the symbol similar. Then saves it
        to proper file.

        Args:
            sample (list of lists of int): list of feature-vectors,
                                           on which we base on.
            symbol (String): name of symbol to compute tolerance
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
        tolerance_distance = means[critical_index] * 1.3
        print("tolerance distance: %.16f" % tolerance_distance)

        tolerance_distance_path = \
            Classifier._get_file_path(
                self.files[DISTANCE_TOLERANCE_FILE], symbol)

        with open(tolerance_distance_path, 'w') as handle:
            handle.write("%.16f\n" % tolerance_distance)

        return tolerance_distance

    def _write_training_set_to_file(self, symbol):
        """Write current training set to file."""

        file_with_training_path = \
            Classifier._get_file_path(self.files[TRAINING_SET_FILE], symbol)

        with open(file_with_training_path, 'wb') as handle:
            pickle.dump(self.training_set, handle)

    def _save_symbol_list(self):
        """Save current list of symbols."""

        with open(self.files[SYMBOL_LIST_FILE], 'wb') as handle:
            pickle.dump(self.symbol_list, handle)

    def _save_training_set(self, symbol):
        """Save the drawn training set to file.

        Args:
            symbol (str): Name of the symbol.
        """

        self._write_training_set_to_file(symbol)

        if symbol not in self.symbol_list:
            self.symbol_list.append(symbol)
            self._save_symbol_list()

    def _get_inactive_symbols(self):
        """Fetch list of inactive symbols from file."""
        try:

            with open(self.files[INACTIVE_SYMBOLS_FILE], 'rb') as handle:
                inactive_symbols = pickle.load(handle)
            print(inactive_symbols)
            return inactive_symbols

        except FileNotFoundError:
            return []

    def _save_inactive_symbols(self, inactive_symbols):
        """Save given inactive symbols list to file.

        Args:
            inactive_symbols (list of str): Current inactive symbols.
        """
        with open(self.files[INACTIVE_SYMBOLS_FILE], 'wb') as handle:
            pickle.dump(inactive_symbols, handle)        

    def _learn_one_symbol(self, symbol):
        """Learn given symbol basing on training set from file.

        Args:
            symbol (str): Name of the symbol.
        """
        training_set = self._load_training_set(symbol)
        feature_vectors = []
        for training_element in training_set:
            feature_vectors.append(featureextractor
                                   .get_features(training_element))
        sample = np.array(feature_vectors)
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree')\
            .fit(sample)

        model_path = Classifier.\
            _get_file_path(self.files[MODEL_FILE], symbol)

        with open(model_path, 'wb') as handle:
            pickle.dump(nbrs, handle)

        return self._compute_tolerance_distance(sample, symbol)

    def _learn_all_symbols_together(self):
        """Build file of knn-classifier model of all training elements."""
        print('learning all together...')
        feature_vectors = []
        results = []
        for sym in self.symbol_list:

            training_set = self._load_training_set(sym)
            for training_element in training_set:
                feature_vector = \
                    featureextractor.get_features(training_element)
                feature_vectors.append(feature_vector)
                results.append(sym)
        if self.symbol_list:
            knn_model = KNeighborsClassifier(n_neighbors=5).\
                fit(feature_vectors, results)
            with open(Classifier._get_file_path(self.files[MODEL_FILE], ""),
                      'wb') as handle:
                pickle.dump(knn_model, handle)

        else:
            try:
                os.remove(Classifier._get_file_path(self.files[MODEL_FILE],
                                                    ""))
            except OSError:
                pass

    def learn(self, load_from_file, symbol=""):
        """Learn basing on training-set.

        Args:
            load_from_file (bool): True - if training has to be load from file,
                          False - new training-set written in self.training_set
                                 that has to be learned and then saved to file.
            symbol (str): Name of the symbol,
                          or empty string if general learning wanted.
        """
        if symbol is None:
            symbol = ""
        if symbol != "":
            print("learning", symbol, "symbol...")
            if not load_from_file:
                self._save_training_set(symbol)
            self._learn_one_symbol(symbol)
        else:
            for sym in self.symbol_list:
                print("learning", sym, "symbol...")
                self._learn_one_symbol(sym)

        self._learn_all_symbols_together()

    def _delete_symbol(self, symbol):
        print('removing symbol', symbol, 'from classifier...')
        inactive_symbols = self._get_inactive_symbols()
        if symbol in self.symbol_list:
            self.symbol_list.remove(symbol)
            with open(self.files[SYMBOL_LIST_FILE], 'wb') as handle:
                handle.truncate()
                pickle.dump(self.symbol_list, handle)

        elif symbol in inactive_symbols:
            inactive_symbols.remove(symbol)
            self._save_inactive_symbols(inactive_symbols)

        else:
            print('warning: symbol', symbol,
                  'is not present in classifier database')
            for symbol in self.symbol_list:
                print(symbol)

        print("removing related files...")
        try:
            os.remove(Classifier._get_file_path(
                self.files[TRAINING_SET_FILE], symbol))
        except OSError:
            pass
        try:
            os.remove(Classifier._get_file_path(
                self.files[MODEL_FILE], symbol))
        except OSError:
            pass
        try:
            os.remove(Classifier._get_file_path(
                self.files[DISTANCE_TOLERANCE_FILE], symbol))
        except OSError:
            pass

    def delete_symbols(self, symbols_to_delete):
        """Delete symbols from classifier with all files related.

        Args:
            symbols_to_delete (list of str): Symbols to delete names.
        """
        if not symbols_to_delete:
            print('removing all symbols from classifier')
            symbols_to_delete = self.symbol_list.copy()

        if symbols_to_delete:
            for symbol in symbols_to_delete:
                self._delete_symbol(symbol)

        self._learn_all_symbols_together()

    def activate_symbols(self, symbols):
        inactive_symbols = self._get_inactive_symbols()
        if not symbols:
            symbols = inactive_symbols.copy()
        allowed = True
        for symbol in symbols:
            if symbol in inactive_symbols:
                print("activating symbol", symbol, "in classifier...")
                if not os.path.isfile(Classifier._get_file_path(self.files[TRAINING_SET_FILE], symbol)):
                    print("File with training set of symbol", symbol, "is missing. Activation is impossible.")
                    allowed = False
                elif not os.path.isfile(Classifier._get_file_path(self.files[MODEL_FILE], symbol)):
                    print("File with learning model of symbol", symbol, "is missing. Activation is impossible.")
                    allowed = False
                elif not os.path.isfile(Classifier._get_file_path(self.files[DISTANCE_TOLERANCE_FILE], symbol)):
                    print("File with the tolerance distance of symbol", symbol, "is missing. Activation is impossible.")
                    allowed = False
                else:
                    self.symbol_list.append(symbol)
                    inactive_symbols.remove(symbol)
            elif symbol not in self.symbol_list:
                print("warning: symbol", symbol, "is not present in classifier database")
        if allowed:
            self._save_symbol_list()
            self._save_inactive_symbols(inactive_symbols)
            self._learn_all_symbols_together()
            print("activation in classifier passed with success")
            return True
        return False

    def deactivate_symbols(self, symbols):
        if not symbols:
            symbols = self.symbol_list.copy()
        inactive_symbols = self._get_inactive_symbols()
        for symbol in symbols:
            if symbol in self.symbol_list:
                print("deactivating symbol", symbol, "in classifier...")
                self.symbol_list.remove(symbol)
                inactive_symbols.append(symbol)
            elif symbol not in inactive_symbols:
                print("warning: symbol", symbol, "is not present in classifier database")
        self._save_symbol_list()
        self._save_inactive_symbols(inactive_symbols)
        self._learn_all_symbols_together()

    @staticmethod
    def _get_file_path(template_string, symbol_name):
        """Transform the file path template to real path.

        Args:
            template_string (str): Template of file path.
            symbol_name (str): Name of symbol or empty string
                               for general files.

        Returns:
            Actual file path.
        """
        subst = "_" + symbol_name
        if symbol_name == "":
            subst = ""

        return Template(template_string).substitute(sym=subst)

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

        file_paths = [operator.add(l, r)
                      for l, r in zip(file_paths, files)]

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
