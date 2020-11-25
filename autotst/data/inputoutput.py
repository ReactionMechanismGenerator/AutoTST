#!/usr/bin/python
# -*- coding: utf-8 -*-

##########################################################################
#
#   AutoTST - Automated Transition State Theory
#
#   Copyright (c) 2015-2020 Richard H. West (r.west@northeastern.edu)
#   and the AutoTST Team
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
from ..reaction import Reaction, TS
from ..species import Species, Conformer
from ..calculator.vibrational_analysis import VibrationalAnalysis
from .base import QMData
import rmgpy.molecule
import rmgpy.species
import rmgpy.reaction
import rmgpy.kinetics
import rmgpy.statmech 


def get_possible_names(reactants, products):

    joined_react_orders = ['+'.join(order)
                         for order in itertools.permutations(reactants)]
    joined_prod_orders = ['+'.join(order)
                        for order in itertools.permutations(products)]
    file_names = [
        '_'.join(
            (r_jo, p_jo)) for r_jo in joined_react_orders for p_jo in joined_prod_orders] + [
        '_'.join(
                (p_jo, r_jo)) for r_jo in joined_react_orders for p_jo in joined_prod_orders]

    return file_names


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
            logging.info(f"Sorry, we cannot find a valid file path for {self.reaction}...")
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
            logging.info(f"No QMData for {self.reaction}...")
            return False
        file_path = self.get_ts_file_path()
        logging.info(f"Saving TS result file {file_path}")
        with open(file_path, 'w') as results_file:
            results_file.write(f'rxn_label = "{self.reaction.label!s}\"\n')
            results_file.write(f'method = "{self.method!s}\"\n')
            results_file.write(f"qm_data = {self.qmdata!r}\n")
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
        file_path = self.get_kinetics_file_path()
        logging.info(f"Saving kinetics data file {file_path}")
        with open(file_path, 'w') as result_file:

            assert self.reaction.rmg_reaction.kinetics, "No kinetics calclated for this reaction..."
            result_file.write(f'method = "{self.method!s}\"\n')
            result_file.write(
                f'reaction = {self.reaction.rmg_reaction!r}\n')
        return True

    def read_ts_file(self, path=None):
        """
        Load the specified transition state data file and return the dictionary of its contents.

        Variables:
        - reaction (Reaction): the reaction of interest

        Returns:
        - local_context (dict): the content of the .ts file. Will return None if there is an error
        """
        if not path:
            path = self.get_ts_file_path()
        if not os.path.exists(path):
            logging.info(f"We could not file a ts file for {self.reaction.label}")
            return None
        try:
            with open(path) as result_file:
                logging.info(f'Reading existing ts file {path}')
                global_context = {'__builtins__': None}
                local_context = {
                    '__builtins__': None,
                    'True': True,
                    'False': False,
                    'QMData': QMData,
                    'array': np.array,
                    'int32': np.int32,
                }
                exec(result_file.read(), global_context, local_context)
        except IOError as e:
            logging.info(f"Couldn't read ts file {path}")
            return None
        except (NameError, TypeError, SyntaxError) as e:
            logging.error(f'The ts file "{path}" was invalid:')
            logging.exception(e)
            return None
        if 'rxn_label' not in local_context:
            logging.error(
                f'The ts file "{path}" did not contain a rxn_label.')
            return None
        if 'method' not in local_context:
            logging.error(
                f'The ts file "{path}" did not contain a method.')
            return None
        if 'qm_data' not in local_context:
            logging.error(
                f'The ts file "{path}" did not contain thermo_data.')
            return None
        return local_context

    def read_kinetics_file(self, path=None):
        """
        Load the specified kinetic data file and return the dictionary of its contents.

        Variables:
        - reaction (Reaction): the reaction of interest

        Returns:
        - local_context (dict): the content of the .kinetics file. Will return None if there is an error
        """
        if not path:
            path = self.get_kinetics_file_path()
        if not os.path.exists(path):
            logging.info(f"We could not file a kinetics file for {self.reaction.label}")
            return None
        try:
            with open(path) as result_file:
                logging.info(f'Reading existing kinetics file {path}')
                global_context = {
                    #'__builtins__': None
                }
                local_context = {
                    #'__builtins__': None,
                    'True': True,
                    'False': False,
                    'Reaction': rmgpy.reaction.Reaction,
                    'TemplateReaction' : rmgpy.data.kinetics.TemplateReaction,
                    'Species': rmgpy.species.Species,
                    'TransitionState': rmgpy.species.TransitionState,
                    'Arrhenius': rmgpy.kinetics.Arrhenius,
                    'Eckart': rmgpy.kinetics.Eckart,
                    'Conformer': rmgpy.statmech.Conformer,
                    'IdealGasTranslation': rmgpy.statmech.IdealGasTranslation,
                    'NonlinearRotor': rmgpy.statmech.NonlinearRotor,
                    'HarmonicOscillator': rmgpy.statmech.HarmonicOscillator,
                    'LinearRotor': rmgpy.statmech.LinearRotor,
                    'array': np.array,
                    'int32': np.int32,
                    'Molecule': rmgpy.molecule.Molecule,
                    'SingleExponentialDown':rmgpy.pdep.SingleExponentialDown
                }
                exec(result_file.read(), global_context, local_context)
        except IOError as e:
            logging.error(f"Couldn't read kinetics file {path}")
            return None
        except (NameError, TypeError, SyntaxError) as e:
            logging.error(f'The kinetics file "{path}" was invalid:')
            logging.exception(e)
            return None
        if 'method' not in local_context:
            logging.error(
                f'The kinetics file "{path}" did not contain a method.')
            return None
        if 'Reaction' not in local_context:
            logging.error(
                f'The kinetics file "{path}" did not contain a reaction.')
            return None
        return local_context
