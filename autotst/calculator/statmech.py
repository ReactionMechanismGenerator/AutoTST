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
import logging
import ase
import cclib.io

from ..reaction import Reaction, TS
from ..species import Species, Conformer

import rmgpy
import arkane.main
from arkane.ess import GaussianLog

FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)


class StatMech():

    def __init__(
            self,
            reaction = None,
            direction = "forward",
            choose_exothermic_direction = True,
            species_directory = "./species",
            ts_directory = "./ts",
            model_chemistry= "M06-2X/cc-pVTZ",
            atom_corrections = False,
            bond_corrections = False,
            hindered_rotors = False,
            freq_scale_factor= 0.982):
        """
        A class to perform Arkane calculations:
        :param: reaction: (Reaction) The reaction of interest
        :param: output_directory: (str) The directory where you would like output files written to
        :param: model_chemistry: (str) The supported model_chemistry described by http://reactionmechanismgenerator.github.io/RMG-Py/users/arkane/input.html#model-chemistry
        :param: freq_scale_factor: (float) The scaling factor corresponding to the model chemistry - source:https://comp.chem.umn.edu/freqscale/version3b1.htm
        """

        if isinstance(reaction,str):
            self.reaction = Reaction(reaction)
        elif isinstance(reaction, Reaction):
            self.reaction = reaction
        else:
            raise TypeError(
                "reaction provided must be class :str: or class :autotst.reaction.Reaction:, not "
                f"{type(reaction)}")

        assert direction.lower() in ("forward","reverse")
        self.direction = direction.lower()
        self.choose_exothermic_direction = choose_exothermic_direction
        self.species_directory = species_directory
        self.ts_directory = ts_directory
        
        self.kinetics_job = arkane.main.Arkane()
        self.thermo_job = arkane.main.Arkane()
        self.model_chemistry = model_chemistry
        self.freq_scale_factor = freq_scale_factor

        self.atom_corrections = atom_corrections
        self.bond_corrections = bond_corrections
        self.hindered_rotors = hindered_rotors

        self.reactants = {}
        self.products = {}
        self.ts = {}

        self.energies = [0.0,0.0,0.0]
        if self.choose_exothermic_direction is True:
            self._calculate_wells()
            self.direction  = self._determine_exothermic_direction()

    def _calculate_e0(self,path):

        if not os.path.exists(path):
            logging.warning(
                f"{path} does not exists...cannot calculate energy")
            return None

        parser = GaussianLog(path)
        e_elect = parser.load_energy()
        try:
            e_zpe = parser.load_zero_point_energy()
        except:
            e_zpe = 0
        E0 = (e_elect + e_zpe) * 1000.0 # kJ/mol

        return E0

    def _calculate_wells(self):


        # Determine energy of reactants
        for react in self.reaction.reactants:
            lowest_energy = 1e5
            for smiles in list(react.conformers.keys()):
                path = os.path.join(self.species_directory,
                                    smiles, smiles + ".log")
                e0 = self._calculate_e0(self, path)
                if e0 is None:
                    continue

                if e0 < lowest_energy:
                    lowest_energy = e0
                    lowest_energy_smiles = smiles
                    lowest_energy_path = path

            self.energies[0] += lowest_energy
            self.reactants[lowest_energy_smiles] = (lowest_energy_path, lowest_energy)

        # Determine energy of TS
        label = self.reaction.label
        path = os.path.join(self.ts_directory, label, label + ".log")
        e0 = self._calculate_e0(self, path)
        self.energies[1] = e0
        self.ts[label] = (path, e0)
        
        # Determine energy of products
        for prod in self.reaction.products:
            lowest_energy = 1e5
            for smiles in list(prod.conformers.keys()):
                path = os.path.join(self.species_directory,
                                    smiles, smiles + ".log")
                e0 = self._calculate_e0(self, path)
                if e0 is None:
                    continue

                if e0 < lowest_energy:
                    lowest_energy = e0
                    lowest_energy_smiles = smiles
                    lowest_energy_path = path

            self.energies[2] += lowest_energy
            self.products[lowest_energy_smiles] = (lowest_energy_path, lowest_energy)

        logging.info(f"The energies of the reactants, ts, and products (in kJ/mol) are {self.energies}")
        e_reacts, e_ts, e_prods = self.energies

        if e_ts <= min(e_reacts,e_prods):
            logging.warning("NEGATIVE BARRIER HEIGHT: "
                            f"The transition state energy ({round(e_ts,3)} kJ/mol) is less than the "
                            f"energies of the reactants ({round(e_reacts,3)} kJ/mol) "
                            f"and the products ({round(e_prods,3)} kJ/mol)")
        
        return self.energies

    def _determine_exothermic_direction(self):

        if all([e == 0.0 for e in self.energies]):
            self._calculate_wells()
        
        e_reacts, e_ts, e_prods = self.energies
        
        if e_reacts > e_prods:
            logging.info(f"The products are lower in energy by {round(e_reacts-e_prods,3)} kJ/mol")
            logging.info("The forward direction is exothermic")
            return "forward"
        else:
            logging.info(f"The reactants are lower in energy by {round(e_prods-e_reacts,3)} kJ/mol")
            logging.info("The reverse direction is exothermic")
            return "reverse"


    def write_species_files(self, species):
        """
        A method to write Arkane files for all conformers in a Species object

        Parameters:
        - species (Species): a species object that you want to write arkane files for
        - scratch (str): the directory where you want to write arkane files to, there should be a 'species/SMILES/' subdirectory

        Returns:
        - None
        """

        for smiles, confs in list(species.conformers.items()):
            if os.path.exists(os.path.join(self.species_directory, smiles, smiles + ".log")):
                logging.info(
                    f"Lowest energy conformer log file exists for {smiles}")
                self.write_conformer_file(conformer=confs[0])
            else:
                logging.info(
                    f"Lowest energy conformer log file DOES NOT exist for {smiles}")

    def write_conformer_file(self, conformer):
        """
        A method to write Arkane files for a single Conformer object

        Parameters:
        - conformer (Conformer): a Conformer object that you want to write an Arkane file for
        - scratch (str): the directory where you want to write arkane files to, there should be a 'species/SMILES/' subdirectory

        Returns:
        - None
        """
        label = conformer.smiles

        if not os.path.exists(os.path.join(self.species_directory, label, label + ".log")):
            logging.info("There is no lowest energy conformer file...")
            return False

        if os.path.exists(os.path.join(self.species_directory, label, label + '.py')):
            logging.info("Species input file already written... Not doing anything")
            return True

        parser = cclib.io.ccread(os.path.join(
            self.species_directory, label, label + ".log"), loglevel=logging.ERROR)
        
        symbol_dict = {
            35: "Br",
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }

        atoms = []

        for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
            atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))

        conformer._ase_molecule = ase.Atoms(atoms)
        conformer.update_coords_from("ase")
        mol = conformer.rmg_molecule
        output = ['#!/usr/bin/env python',
                  '# -*- coding: utf-8 -*-', ]

        output += ["",
                   f"spinMultiplicity = {parser.mult}",
                   ""]

        output += ["energy = {", f"    '{self.model_chemistry}': Log('{label}.log'),", "}", ""]  # fix this

        output += [f"geometry = Log('{label}.log')", ""]

        output += [
            f"frequencies = Log('{label}.log')", ""]

        """
        TODO: add rotor information @carl
        output += ["rotors = ["]
        for torsion in conf.torsions:
            output += [self.get_rotor_info(conf, torsion)]
        output += ["]"]
        """
        input_string = ""

        for t in output:
            input_string += t + "\n"

        with open(os.path.join(self.species_directory, label, label + '.py'), "w") as f:
            f.write(input_string)
        return True

    def get_rotor_info(self, conformer, torsion_index):
        """
        Formats and returns info about torsion as it should appear in an Arkane species.py

        The following are needed for an Arkane input file:
        - scanLog :: Gaussian output log of freq calculation on optimized geometry
        - pivots  :: torsion center: j,k of i,j,k,l (Note Arkane begins indexing with 1)
        - top     :: ID of all atoms in one top (Note Arkane begins indexing with 1)

        Parameters:
        - conformer (Conformer): autotst conformer object
        - torsion (Torsion): autotst torsion object 


        Returns:
        - info (str): a string containing all of the relevant information for a hindered rotor scan
        """
        torsion = conformer.torsions[torsion_index]
        _, j, k, _ = torsion.atom_indices

        # Adjusted since mol's IDs start from 0 while Arkane's start from 1
        tor_center_adj = [j+1, k+1]

        if isinstance(conformer, TS):
            tor_log = os.path.join(
                scratch,
                "ts",
                conformer.reaction_label,
                "torsions",
                comformer.reaction_label + f"_36by10_{j}_{k}.log"
            )
        else:
            tor_log = os.path.join(
                scratch,
                "species",
                conformer.smiles,
                "torsions",
                conformer.smiles + f'_36by10_{j}_{k}.log'
            )

        if not os.path.exists(tor_log):
            logging.info(
                f"Torsion log file does not exist for {torsion}")
            return ""

        top_IDs = []
        for num, tf in enumerate(torsion.mask):
            if tf:
                top_IDs.append(num)

        # Adjusted to start from 1 instead of 0
        top_IDs_adj = [ID+1 for ID in top_IDs]

        info = f"     HinderedRotor(scanLog=Log('{tor_log}'), pivots={tor_center_adj}, top={top_IDs_adj}, fit='fourier'),"

        return info

    def write_ts_input(self):
        """
        A method to write Arkane files for a single TS object

        Parameters:
        - transitionstate (TS): a TS object that you want to write an Arkane file for
        - scratch (str): the directory where you want to write arkane files to, there should be a 'ts/REACTION_LABEL/' subdirectory

        Returns:
        - None
        """

        label = self.reaction.label

        if os.path.exists(os.path.join(self.ts_directory, label, label + '.py')):
            logging.info("TS input file already written... Not doing anything")
            return True

        if not os.path.exists(os.path.join(self.ts_directory, label, label + ".log")):
            logging.info("There is no lowest energy conformer file...")
            return False

        parser = cclib.io.ccread(os.path.join(self.ts_directory, label, label + ".log"), loglevel=logging.ERROR)
        
        symbol_dict = {
            35: "Br",
            17: "Cl",
            9:  "F",
            8:  "O",
            7:  "N",
            6:  "C",
            1:  "H",
        }


        # atoms = []
        # for atom_num, coords in zip(parser.atomnos, parser.atomcoords[-1]):
        #     atoms.append(ase.Atom(symbol=symbol_dict[atom_num], position=coords))

        # transitionstate = self.reaction.ts["forward"][0]
        # transitionstate._ase_molecule = ase.Atoms(atoms)
        # transitionstate.update_coords_from("ase")

        output = ['#!/usr/bin/env python',
                  '# -*- coding: utf-8 -*-']

        # transitionstate.rmg_molecule.update_multiplicity()

        output += ["",
                   f"spinMultiplicity = {parser.mult}",
                   ""]

        output += ["energy = {", f"    '{self.model_chemistry}': Log('{label}.log'),", "}", ""]  # fix this

        output += [f"geometry = Log('{label}.log')", ""]

        output += [
            f"frequencies = Log('{label}.log')", ""]
        
        if self.hindered_rotors is True:
            output += ["rotors = []", ""]  # TODO: Fix this

        input_string = ""

        for t in output:
            input_string += t + "\n"

        with open(os.path.join(self.ts_directory, label, label + '.py'), "w") as f:
            f.write(input_string)
        return True

    def write_kinetics_input(self):
        """
        A method to write Arkane file to obtain kinetics for a Reaction object

        Parameters:
        - reaction (Reaction): a Reaction object that you want to write an Arkane kinetics job file for.
        - scratch (str): the directory where you want to write arkane files to, there should be a 'species/SMILES/' subdirectory

        Returns:
        - None
        """

        top = [
            "#!/usr/bin/env python",
            "# -*- coding: utf-8 -*-",
            "",
            f'modelChemistry = "{self.model_chemistry}"',  # fix this
            f"frequencyScaleFactor = {self.freq_scale_factor}",  # fix this
            f"useAtomCorrections = {self.atom_corrections}"
            f"useHinderedRotors = {self.hindered_rotors}",  # fix this @carl
            f"useBondCorrections = {self.bond_corrections}",
            ""]

        if len(self.reactants) == 0:
            self._calculate_wells()
        
        r_labels = []
        p_lables = []
        if self.direction == "forward":
            for i, (smiles, (path, energy)) in enumerate(self.reactants.items()):
                label = f"react_{i}"
                line = f"species('{label}', '{path}', structure=SMILES('{smiles}'))"
                top.append(line)
                r_labels.append(label)
            for i, (smiles, (path, energy)) in enumerate(self.products.items()):
                label = f"prod_{i}"
                line = f"species('{label}', '{path}', structure=SMILES('{smiles}'))"
                top.append(line)
                p_labels.append(label)
        else: # direction == 'reverse'
            for i, (smiles, (path, energy)) in enumerate(self.products.items()):
                label = f"react_{i}"
                line = f"species('{label}', '{path}', structure=SMILES('{smiles}'))"
                top.append(line)
                r_labels.append(label)
            for i, (smiles, (path, energy)) in enumerate(self.reactants.items()):
                label = f"prod_{i}"
                line = f"species('{label}', '{path}', structure=SMILES('{smiles}'))"
                top.append(line)
                p_labels.append(label)
    
        ts_path = os.path.join(self.ts_directory, self.reaction.label, self.reaction.label + ".py")
        line = f"transitionState('TS', '{ts_path}')"
        top.append(line)

        line = ["",
                "reaction(",
                f"    label = '{self.reaction.label}',",
                f"    reactants = {r_labels},",
                f"    products = {p_labels},",
                "    transitionState = 'TS',",
                "    tunneling = 'Eckart',",
                ")",
                "",
                "statmech('TS')",
                f"kinetics('{self.reaction.label}')"]

        top += line

        input_string = ""

        for t in top:
            input_string += t + "\n"

        with open(os.path.join(self.ts_directory, self.reaction.label, self.reaction.label + ".kinetics.py"), "w") as f:
            f.write(input_string)

    def write_thermo_input(self, conformer):
        """
        A method to write Arkane file to obtain thermochemistry for a Conformer object

        Parameters:
        - conformer (Conformer): a Conformer object that you want to write an Arkane thermo job file for.
        - scratch (str): the directory where you want to write arkane files to, there should be a 'species/SMILES/' subdirectory

        Returns:
        - None
        """

        model_chemistry = self.model_chemistry
        freq_scale_factor = self.freq_scale_factor

        top = [
            "#!/usr/bin/env python",
            "# -*- coding: utf-8 -*-",
            "",
            f'modelChemistry = "{model_chemistry}"',  # fix this
            f"frequencyScaleFactor = {freq_scale_factor}",  # fix this
            "useHinderedRotors = False",  # fix this @carl
            "useBondCorrections = False",
            ""]
        p = os.path.join(conformer.smiles + ".py")
        line = f"species('species', '{f}', structure=SMILES('{conformer.smiles}'))"
        top.append(line)

        top.append("statmech('species')")
        top.append("thermo('species', 'NASA')")

        input_string = ""

        for t in top:
            input_string += t + "\n"

        with open(os.path.join(self.species_directory, conformer.smiles, conformer.smiles + ".thermo.py"), "w") as f:
            f.write(input_string)

    def write_files(self):
        """
        A method to write all species, transition state, and kinetics job files to obtain kinetic parameters

        Parameters:
        - None

        Returns:
        - None
        """

        for mol in self.reaction.reactants + self.reaction.products:
            for smiles, confs in mol.conformers.items:
                conf = confs[0]
                self.species[smiles] = conf
                self.write_conformer_file(conf)

        self.write_ts_input(
            self.reaction.ts["forward"][0])

        self.write_kinetics_input()

    def run(self):
        """
        A method to write run a kinetics job from all of the files written `write_files`

        Parameters:
        - None

        Returns:
        - None
        """

        self.kinetics_job.input_file = os.path.join(
            self.ts_directory, self.reaction.label, self.reaction.label + ".kinetics.py")
        self.kinetics_job.plot = False
        self.kinetics_job.output_directory = os.path.join(self.ts_directory, self.reaction.label)

        self.kinetics_job.execute()

        for job in self.kinetics_job.job_list:
            if isinstance(job, arkane.main.KineticsJob):
                self.kinetics_job = job
            elif isinstance(job, arkane.main.ThermoJob):
                self.thermo_job = job

    def set_results(self):
        """
        A method to set the RMGReaction from the kinetics job to the RMGReaction of the input Reaction

        Parameters:
        - None

        Returns:
        - None
        """

        for reactant in self.reaction.rmg_reaction.reactants:
            for r in self.kinetics_job.reaction.reactants:
                if reactant.to_smiles() == r.label:
                    r.molecule = [reactant]

        for product in self.reaction.rmg_reaction.products:
            for p in self.kinetics_job.reaction.products:
                if product.to_smiles() == p.label:
                    p.molecule = [product]

        self.reaction.rmg_reaction = self.kinetics_job.reaction

        return self.reaction
