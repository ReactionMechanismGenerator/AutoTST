#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2018 Prof. Richard H. West (r.west@northeastern.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
##########################################################################

import os
import itertools
import logging
import numpy as np
from autotst.reaction import Reaction, TS
from autotst.species import Species, Conformer
from autotst.calculator.vibrational_analysis import VibrationalAnalysis
from autotst.data.base import QMData
import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.molecule import Atom, getElement
from rmgpy.species import Species as RMGSpecies, TransitionState
from rmgpy.reaction import Reaction as RMGReaction
from rmgpy.kinetics import Arrhenius, Eckart
from rmgpy.statmech import Conformer as RMGConformer, IdealGasTranslation, NonlinearRotor, HarmonicOscillator, LinearRotor


def get_possible_names(reactants, products):

    joinedReactOrders = ['+'.join(order)
                         for order in itertools.permutations(reactants)]
    joinedProdOrders = ['+'.join(order)
                        for order in itertools.permutations(products)]
    fileNames = [
        '_'.join(
            (rJO, pJO)) for rJO in joinedReactOrders for pJO in joinedProdOrders] + [
        '_'.join(
                (pJO, rJO)) for rJO in joinedReactOrders for pJO in joinedProdOrders]

    return fileNames


class InputOutput():
    """
    A class that allows users to read and write .ts and .kinetics files for validated TSs and Reactions
    """

    def __init__(self, reaction, method="m062x", directory="."):
        """
        Variables:
        - reaction (Reaction): the reaction of interest
        - directory (str): the directory where you want .ts and .kinetics files saved
        """
        self.method = method
        self.reaction = reaction
        self.label = self.reaction.label
        self.directory = directory

    def copy(self):
        """
        A method to generate a copy of this object
        """
        from copy import deepcopy
        return deepcopy(self)

    def get_qmdata(self):
        """
        A method to create a QMData object from a log file at file_path

        Variables:
        - file_path (str): the location of the log file you want in
        """

        file_path = os.path.join(self.directory, "ts", self.reaction.label, self.reaction.label + ".log")
        if not os.path.exists(file_path):
            logging.info("Sorry, we cannot find a valid file path for {}...".format(self.reaction))
            self.qmdata = None
            return self.qmdata
        self.qmdata = QMData()
        self.qmdata.get_qmdata(file_path=file_path)
        return self.qmdata

    def get_ts_file_path(self):

        return os.path.join(self.directory, "ts", self.reaction.label, self.reaction.label + ".ts")

    def get_kinetics_file_path(self):

        return os.path.join(self.directory, "ts", self.reaction.label, self.reaction.label + ".kinetics")

    def save_ts(self):
        """
        A method to save the .ts data to a file

        Variables:
        - reaction (Reaction): the reaction of interest
        - method (str): the computational methods that were used
        - qmData (QMData): the data that you want to save
        """
        self.get_qmdata()
        if not self.qmdata:
            logging.info("No QMData for {}...".format(self.reaction))
            return False
        file_path = self.get_ts_file_path()
        logging.info("Saving TS result file {}".format(file_path))
        with open(file_path, 'w') as resultFile:
            resultFile.write('rxnLabel = "{0!s}"\n'.format(self.reaction.label))
            resultFile.write('method = "{0!s}"\n'.format(self.method))
            resultFile.write("qmData = {0!r}\n".format(self.qmdata))
        return True

    def save_kinetics(self):
        """
        A method to save the .ts data to a file

        Variables:
        - reaction (Reaction): the reaction of interest
        - method (str): the computational methods that were used
        """
        if not self.reaction.rmg_reaction.kinetics:
            logging.info("The reaction you're trying to save info for doesn't have any kinetic information... aborting.")
            return False
        filePath = self.get_kinetics_file_path()
        logging.info("Saving kinetics data file {}".format(filePath))
        with open(filePath, 'w') as resultFile:

            assert self.reaction.rmg_reaction.kinetics, "No kinetics calclated for this reaction..."
            resultFile.write('method = "{0!s}"\n'.format(self.method))
            resultFile.write(
                'reaction = {0!r}\n'.format(self.reaction.rmg_reaction))
        return True

    def read_ts_file(self):
        """
        Load the specified transition state data file and return the dictionary of its contents.

        Variables:
        - reaction (Reaction): the reaction of interest

        Returns:
        - local_context (dict): the content of the .ts file. Will return None if there is an error
        """

        ts_file = self.get_ts_file_path()
        if os.path.exists(ts_file):
            path = ts_file
        else:
            logging.info("We could not file a ts file for {}".format(self.reaction.label))
            return None
        try:
            with open(path) as resultFile:
                logging.info('Reading existing ts file {0}'.format(path))
                global_context = {'__builtins__': None}
                local_context = {
                    '__builtins__': None,
                    'True': True,
                    'False': False,
                    'QMData': QMData,
                    'array': np.array,
                    'int32': np.int32,
                }
                exec resultFile in global_context, local_context
        except IOError as e:
            logging.info("Couldn't read ts file {0}".format(path))
            return None
        except (NameError, TypeError, SyntaxError) as e:
            logging.error('The ts file "{0}" was invalid:'.format(path))
            logging.exception(e)
            return None
        if 'rxnLabel' not in local_context:
            logging.error(
                'The ts file "{0}" did not contain a rxnLabel.'.format(path))
            return None
        if 'method' not in local_context:
            logging.error(
                'The ts file "{0}" did not contain a method.'.format(path))
            return None
        if 'qmData' not in local_context:
            logging.error(
                'The ts file "{0}" did not contain thermoData.'.format(path))
            return None
        return local_context

    def read_kinetics_file(self):
        """
        Load the specified kinetic data file and return the dictionary of its contents.

        Variables:
        - reaction (Reaction): the reaction of interest

        Returns:
        - local_context (dict): the content of the .kinetics file. Will return None if there is an error
        """

        ts_file = self.get_kinetics_file_path()
        if os.path.exists(ts_file):
            path = ts_file
        else:
            logging.info("We could not file a kinetics file for {}".format(self.reaction.label))
            return None
        try:
            with open(path) as resultFile:
                logging.info('Reading existing kinetics file {0}'.format(path))
                global_context = {'__builtins__': None}
                local_context = {
                    '__builtins__': None,
                    'True': True,
                    'False': False,
                    'Reaction': RMGReaction,
                    'Species': RMGSpecies,
                    'TransitionState': TransitionState,
                    'Arrhenius': Arrhenius,
                    'Eckart': Eckart,
                    'Conformer': RMGConformer,
                    'IdealGasTranslation': IdealGasTranslation,
                    'NonlinearRotor': NonlinearRotor,
                    'HarmonicOscillator': HarmonicOscillator,
                    'LinearRotor': LinearRotor,
                    'array': np.array,
                    'int32': np.int32,
                    'Molecule': RMGMolecule
                }
                exec resultFile in global_context, local_context
        except IOError as e:
            logging.error("Couldn't read kinetics file {0}".format(path))
            return None
        except (NameError, TypeError, SyntaxError) as e:
            logging.error('The kinetics file "{0}" was invalid:'.format(path))
            logging.exception(e)
            return None
        if 'method' not in local_context:
            logging.error(
                'The kinetics file "{0}" did not contain a method.'.format(path))
            return None
        if 'Reaction' not in local_context:
            logging.error(
                'The kinetics file "{0}" did not contain a reaction.'.format(path))
            return None
        return local_context
