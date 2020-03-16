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
import cclib.io

import autotst
from ..reaction import Reaction, TS
from ..species import Species, Conformer
from ..geometry import Torsion

import ase
import ase.calculators.gaussian

import rmgpy
import rmgpy.molecule
import rmgpy.reaction 

class Gaussian():

    def __init__(self,
                 conformer=None, # Either a transition state or a conformer
                 settings={
                        "method": "m062x",
                        "basis": "cc-pVTZ",
                        "dispersion": "",
                        "mem": "5GB",
                        "sp": "",
                        "convergence": "",
                        "nprocshared": 24,
                        "time": "24:00:00",
                 },
                 directory=".", #where you want input and log files to be written, default is current directory
                 scratch=None  #where you want temporary files to be written
                 ):

        default_settings = {
            "method": "m062x",
            "basis": "cc-pVTZ",
            "dispersion": "",
            "mem": "5GB",
            "sp": "",
            "convergence": "",
            "nprocshared": 24,
            "time": "24:00:00",

        }

        self.conformer = conformer

        # setting the settings accordingly
        for setting, value in list(default_settings.items()):
            if setting in list(settings.keys()):
                assert isinstance(settings[setting], type(
                    value)), "{} is not a proper instance..."
            else:
                logging.info("{} not specified, setting it to the default value of {}".format(
                    setting, value))
                settings[setting] = value

        self.command = "g16"
        self.settings = settings

        # Assert convergence is a valid option
        convergence_options = ["", "verytight", "tight", "loose"]
        convergence = self.settings["convergence"].lower()
        assert convergence in convergence_options, f"{convergence} is not among the supported convergence options {convergence_options}"
        
        if settings["dispersion"]:
            dispersion = settings["dispersion"].upper()
            assert dispersion in ['GD3','GD3BJ','GD2'],'Acceptable keywords for dispersion are GD3, GD3BJ, or GD2'
            self.opt_method = settings["method"].upper() + '-' + dispersion + '_' + settings["basis"].upper()
        else:
            self.opt_method = settings["method"].upper() + '_' + settings["basis"].upper()
        
        # Assert disperion is a valid option
        dispersion_options = ['GD3', 'GD3BJ', 'GD2']
        dispersion = self.settings["dispersion"].upper()
        if dispersion is not "":
            assert dispersion in dispersion_options, f"{dispersion} is not among the supported convergence options {dispersion_options}"
            self.settings["dispersion"] = f'EmpiricalDispersion={dispersion}'
            self.model_chemistry = settings["method"].upper() + '-' + dispersion + '_' + settings["basis"].upper()
        else:
            self.model_chemistry = settings["method"].upper() + '_' + settings["basis"].upper()

        self.directory = directory

        try: 
            if scratch is None:
                self.scratch = os.environ['GAUSS_SCRDIR']
            else:
                self.scratch = os.environ['GAUSS_SCRDIR'] = scratch
        except:
            if scratch is None:
                self.scratch = '.'
            else:
                self.scratch = os.environ['GAUSS_SCRDIR'] = scratch

    def __repr__(self):
        if isinstance(self.conformer, TS):
            return f'<Gaussian Calculator {self.conformer.reaction_label}>'
        elif isinstance(self.conformer, Conformer):
            return f'<Gaussian Calculator {self.conformer.smiles}>'
        else:
            return '<Gaussian Calculator>'

    def get_rotor_calc(self,
                       torsion_index=0,
                       steps=36,
                       step_size=10.0):
        """
        A method to create all of the calculators needed to perform hindered rotor calculations given a `Conformer` and a `Torsion`.

        Parameters:
        - conformer (Conformer): A `Conformer` object that you want to perform hindered rotor calculations on
        - torsion (Torsion): A `Torsion` object that you want to perform hindered rotor calculations about
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - scratch (str): a directory where you want log files to be written to
        - steps (int): the number of steps you want performed in this scan
        - step_size (float): the size, in degrees, of the step you to scan along

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """
        torsion = self.conformer.torsions[torsion_index]

        assert (torsion and (isinstance(torsion, Torsion))
                ), "To create a rotor calculator, you must provide a Torsion object."

        assert isinstance(
            self.conformer, Conformer), "A Conformer object was not provided..."

        # Determine number of heavy atoms
        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        time = "24:00:00"

        if num_atoms <= 2:
            mem = "10GB"
            nprocshared = 2
        elif num_atoms <= 4:
            mem = "20GB"
            nprocshared = 6
        elif num_atoms <= 10:
            mem = '40GB'
            nprocshared = 12
        elif num_atoms <= 20:
            mem = '50GB'
            nprocshared = 16
        else:
            mem = '60GB'
            nprocshared = 20

        addsec = ""
        for bond in self.conformer.bonds:
            i, j = bond.atom_indices
            addsec += f"B {(i + 1)} {(j + 1)}\n"

        i, j, k, l = torsion.atom_indices
        addsec += "D {} {} {} {} S {} {}\n".format(
            i + 1, j + 1, k + 1, l + 1, steps, float(step_size))

        if isinstance(self.conformer, TS):
            extra = f"{self.settings['dispersion']} Opt=(CalcFC,ModRedun,{self.settings['convergence']})"
            label = conformer_dir = self.conformer.reaction_label
            label += f"_{steps}by{int(step_size)}_{j}_{k}"
            conformer_type = "ts"

            # Freezing the reaction center for TS HR scans, this is how Pierre did it when running
            # benchmark calcualations for the original AutoTST paper @nateharms
            ind1 = self.conformer.rmg_molecule.get_labeled_atoms("*1")[0].sorting_label
            ind2 = self.conformer.rmg_molecule.get_labeled_atoms("*2")[0].sorting_label
            ind3 = self.conformer.rmg_molecule.get_labeled_atoms("*3")[0].sorting_label
            addsec += f"{(ind1 + 1)} {(ind2 + 1)} F\n"
            addsec += f"{(ind2 + 1)} {(ind3 + 1)} F\n"
            addsec += f"{(ind1 + 1)} {(ind2 + 1)} {(ind3 + 1)} F"

        elif isinstance(self.conformer, Conformer):
            label = conformer_dir = self.conformer.smiles
            label += f"_{steps}by{int(step_size)}_{j}_{k}"
            conformer_type = "species"
            extra = f"{self.settings['dispersion']} Opt=(CalcFC,ModRedun,{self.settings['convergence']})"

        # I don't think we should be freezing the other torsions during a scan
        # for locked_torsion in self.conformer.torsions:  # TODO: maybe doesn't work;
        #     if sorted(locked_torsion.atom_indices) != sorted(torsion.atom_indices):
        #         a, b, c, d = locked_torsion.atom_indices
        #         addsec += f'D {(a + 1)} {(b + 1)} {(c + 1)} {(d + 1)} F\n'

        # Dont think this is needed
        #self.conformer.rmg_molecule.update_multiplicity()
        
        mult = self.conformer.rmg_molecule.multiplicity

        directory = os.path.join(
            self.directory,
            conformer_type,
            conformer_dir,
            "rotors"
        )

        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra=extra,
            multiplicity=mult,
            addsec=[addsec])
        
        ase_gaussian.directory = directory
        ase_gaussian.parameters["partition"] = "" # Create parition parameter to assign in Job
        ase_gaussian.parameters["time"] = time
        ase_gaussian.atoms = self.conformer.ase_molecule
        del ase_gaussian.parameters['force']
        return ase_gaussian

    def get_conformer_calc(self, opt=True, freq=True):
        """
        A method that creates a calculator for a `Conformer` that will perform a geometry optimization

        Parameters:
        - conformer (Conformer): A `Conformer` object that you want to perform hindered rotor calculations on
        - torsion (Torsion): A `Torsion` object that you want to perform hindered rotor calculations about
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - scratch (str): a directory where you want log files to be written to
        - convergence (str): ['verytight','tight','' (default)], specifies the convergence criteria of the geometry optimization

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """


        if isinstance(self.conformer, TS):
            logging.info(
                "TS object provided, cannot obtain a species calculator for a TS")
            return None

        assert isinstance(
            self.conformer, Conformer), "A Conformer object was not provided..."

        #self.conformer.rmg_molecule.update_multiplicity()

        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        if num_atoms <= 4:
            mem = '20GB'
            nprocshared = 16
            time = '1:00:00'
        elif num_atoms <= 6:
            mem = '30GB'
            nprocshared = 24
            time = '1:00:00'
        elif num_atoms <= 12:
            mem = '40GB'
            nprocshared = 8
            time = '12:00:00'
        elif num_atoms <= 20:
            mem = '50GB'
            nprocshared = 12
            time = '12:00:00'
        else:
            mem = '50GB'
            nprocshared = 16
            time = '12:00:00'

        if opt is False:
            time = '1:00:00'

        label = f"{self.conformer.smiles}_{self.conformer.index}"

        new_scratch = os.path.join(
            self.directory,
            "species",
            self.conformer.smiles,
            "conformers"
        )

        try:
            os.makedirs(new_scratch)
        except OSError:
            pass

        extra = f"{self.settings['dispersion']} IOP(7/33=1,2/16=3) scf=(maxcycle=900)"

        if opt is True:
            extra += f" opt=(calcfc,maxcycles=900,{self.settings['convergence']})"
        if freq is True:
            extra += " freq"


        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra=extra,
            multiplicity=self.conformer.rmg_molecule.multiplicity)

        ase_gaussian.directory = new_scratch
        ase_gaussian.parameters["partition"] = ""
        ase_gaussian.parameters["time"] = time
        ase_gaussian.atoms = self.conformer.ase_molecule
        del ase_gaussian.parameters['force']
        return ase_gaussian

    def get_sp_calc(self):

        assert isinstance(
            self.conformer, Conformer), "A Conformer object was not provided..."

        method = self.settings["sp"].upper()

        gaussian_methods = [
            "G1", "G2", "G3", "G4", "G2MP2", "G3MP2", "G3B3", "G3MP2B3", "G4", "G4MP2",
            "W1", "W1U", "W1BD", "W1RO",
            "CBS-4M", "CBS-QB3", "CBS-APNO",
        ]
        assert method in gaussian_methods

        time = "24:00:00"
        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        if num_atoms <= 6:
            mem = '120GB'
            nprocshared = 12
        elif num_atoms <= 8:
            mem = '180GB'
            nprocshared = 16
        elif num_atoms <= 10:
            mem = '250GB'
            nprocshared = 28
        else:
            mem = '380GB'
            nprocshared = 16

        label = "{}_{}".format(self.conformer.smiles, method)

        new_scratch = os.path.join(
            self.directory,
            "species",
            self.conformer.smiles,
            "sp"
        )

        try:
            os.makedirs(new_scratch)
        except OSError:
            pass

        ase_gaussian = ASEGaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=method,
            basis='',
            extra=f"opt=(calcfc,maxcycles=900,{self.settings['convergence']}) IOP(7/33=1,2/16=3) scf=(maxcycle=900)",
            multiplicity=self.conformer.rmg_molecule.multiplicity)
        
        ase_gaussian.atoms = self.conformer.ase_molecule
        ase_gaussian.directory = new_scratch
        ase_gaussian.label = label
        ase_gaussian.parameters["partition"] = ""
        ase_gaussian.parameters["time"] = time
        del ase_gaussian.parameters['force']
        return ase_gaussian

    def get_shell_calc(self):
        """
        A method to create a calculator that optimizes the reaction shell of a `TS` object

        Parameters:
        - ts (TS): A `TS` object that you want to perform calculations on
        - direction (str): the forward or reverse direction of the `TS` object
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - directory (str): a directory where you want log files to be written to

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """
        assert isinstance(self.conformer, TS), "A TS object was not provided..."
        assert self.conformer.direction.lower() in ["forward", "reverse"]

        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        if num_atoms <= 4:
            mem = '20GB'
            nprocshared = 16
            time = '1:00:00'
        elif num_atoms <= 6:
            mem = '30GB'
            nprocshared = 24
            time = '1:00:00'
        elif num_atoms <= 12:
            mem = '40GB'
            nprocshared = 8
            time = '12:00:00'
        elif num_atoms <= 20:
            mem = '50GB'
            nprocshared = 12
            time = '12:00:00'
        else:
            mem = '50GB'
            nprocshared = 16
            time = '12:00:00'

        self.conformer.rmg_molecule.update_multiplicity()

        label = self.conformer.reaction_label + "_" + self.conformer.direction.lower() + "_shell_" + str(self.conformer.index)

        new_scratch = os.path.join(
            self.directory,
            "ts",
            self.conformer.reaction_label,
            "conformers"
        )

        try:
            os.makedirs(new_scratch)
        except OSError:
            pass
        
        #if self.conformer.reaction_family != "Some reaction family with 4 labeled atoms..."
        if self.conformer.reaction_family.lower() in ["h_abstraction", "intra_h_migration", "r_addition_multiplebond"]:
            ind1 = self.conformer.rmg_molecule.get_labeled_atoms("*1")[0].sorting_label
            ind2 = self.conformer.rmg_molecule.get_labeled_atoms("*2")[0].sorting_label
            ind3 = self.conformer.rmg_molecule.get_labeled_atoms("*3")[0].sorting_label
        else:
            logging.error(f"Reaction family {self.conformer.reaction_family} is not supported...")
            raise AssertionError

        combos = ""
        combos += f"{(ind1 + 1)} {(ind2 + 1)} F\n"
        combos += f"{(ind2 + 1)} {(ind3 + 1)} F\n"
        combos += f"{(ind1 + 1)} {(ind2 + 1)} {(ind3 + 1)} F"

        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra=f"{self.settings['dispersion']} Opt=(ModRedun,Loose,maxcycles=900, {self.settings['convergence']}) Int(Grid=SG1) scf=(maxcycle=900)",
            multiplicity=self.conformer.rmg_molecule.multiplicity,
            addsec=[combos]
        )

        ase_gaussian.directory = new_scratch
        ase_gaussian.parameters["time"] = time
        ase_gaussian.parameters["partition"] = ""
        ase_gaussian.atoms = self.conformer.ase_molecule
        del ase_gaussian.parameters['force']

        return ase_gaussian

    def get_center_calc(self):
        """
        A method to create a calculator that optimizes the reaction center of a `TS` object

        Parameters:
        - ts (TS): A `TS` object that you want to perform calculations on
        - direction (str): the forward or reverse direction of the `TS` object
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - scratch (str): a directory where you want log files to be written to

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """

        assert self.conformer.direction.lower() in ["forward", "reverse"]

        assert isinstance(self.conformer, TS), "A TS object was not provided..."

        indicies = []
        for i, atom in enumerate(self.conformer.rmg_molecule.atoms):
            if not (atom.label != ""):
                indicies.append(i)

        addsec = ""
        for combo in list(itertools.combinations(indicies, 2)):
            a, b = combo
            addsec += f"{(a + 1)} {(b + 1)} F\n"

        self.conformer.rmg_molecule.update_multiplicity()

        label = self.conformer.reaction_label + "_" + self.conformer.direction.lower() + "_center_" + str(self.conformer.index)

        new_scratch = os.path.join(
            self.directory,
            "ts",
            self.conformer.reaction_label,
            "conformers"
        )

        try:
            os.makedirs(new_scratch)
        except OSError:
            pass

        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=self.settings["mem"],
            nprocshared=self.settings["nprocshared"],
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra="Opt=(ts,calcfc,noeigentest,ModRedun,maxcycles=900) scf=(maxcycle=900)",
            multiplicity=self.conformer.rmg_molecule.multiplicity,
            addsec=[addsec[:-1]]
        )
        ase_gaussian.atoms = self.conformer.ase_molecule
        del ase_gaussian.parameters['force']

        return ase_gaussian

    def get_overall_calc(self):
        """
        A method to create a calculator that optimizes the overall geometry of a `TS` object

        Parameters:
        - ts (TS): A `TS` object that you want to perform calculations on
        - direction (str): the forward or reverse direction of the `TS` object
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - scratch (str): a directory where you want log files to be written to

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """

        assert isinstance(self.conformer, TS), "A TS object was not provided..."

        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        if num_atoms <= 4:
            mem = '20GB'
            nprocshared = 24
            time = '1:00:00'
        elif num_atoms <= 12:
            mem = '40GB'
            nprocshared = 8
            time = '24:00:00'
        elif num_atoms <= 20:
            mem = '50GB'
            nprocshared = 12
            time = '24:00:00'
        else:
            mem = '50GB'
            nprocshared = 16
            time = '24:00:00'

        self.conformer.rmg_molecule.update_multiplicity()

        label = self.conformer.reaction_label + "_" + self.conformer.direction.lower() + "_" + str(self.conformer.index)

        new_scratch = os.path.join(
            self.directory,
            "ts",
            self.conformer.reaction_label,
            "conformers"
        )

        try:
            os.makedirs(new_scratch)
        except OSError:
            pass

        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra=f"{self.settings['dispersion']} opt=(ts,calcfc,noeigentest,maxcycles=900) freq scf=(maxcycle=900) IOP(7/33=1,2/16=3)",
            multiplicity=self.conformer.rmg_molecule.multiplicity)

        ase_gaussian.directory = new_scratch
        ase_gaussian.parameters["time"] = time
        ase_gaussian.parameters["partition"] = ""
        ase_gaussian.atoms = self.conformer.ase_molecule
        del ase_gaussian.parameters['force']

        return ase_gaussian

    def get_irc_calc(self):
        """
        A method to create a calculator that runs an IRC calculation the overall geometry of a `TS` object

        Parameters:
        - ts (TS): A `TS` object that you want to perform calculations on
        - direction (str): the forward or reverse direction of the `TS` object
        - settings (dict): a dictionary of settings containing method, basis, mem, nprocshared
        - scratch (str): a directory where you want log files to be written to

        Returns:
        - calc (ase.calculators.gaussian.Gaussian): an ase.calculators.gaussian.Gaussian calculator with all of the proper setting specified
        """

        assert isinstance(self.conformer, TS), "A TS object was not provided..."

        num_atoms = self.conformer.rmg_molecule.get_num_atoms(
        ) - self.conformer.rmg_molecule.get_num_atoms('H')

        time = '24:00:00'

        if num_atoms <= 4:
            mem = '20GB'
            nprocshared = 8
        elif num_atoms <= 12:
            mem = '40GB'
            nprocshared = 16
        elif num_atoms <= 20:
            mem = '50GB'
            nprocshared = 20
        else:
            mem = '60GB'
            nprocshared = 24

        self.conformer.rmg_molecule.update_multiplicity()
        label = self.conformer.reaction_label + "_irc_" + self.conformer.direction + "_" + str(self.conformer.index)

        new_scratch = os.path.join(
            self.directory,
            "ts",
            self.conformer.reaction_label,
            "irc"
        )
        try:
            os.makedirs(new_scratch)
        except OSError:
            pass

        ase_gaussian = ase.calculators.gaussian.Gaussian(
            mem=mem,
            nprocshared=nprocshared,
            label=label,
            scratch=new_scratch,
            method=self.settings["method"],
            basis=self.settings["basis"],
            extra=f"{self.settings['dispersion']} irc=(calcall)",
            multiplicity=self.conformer.rmg_molecule.multiplicity
        )
        
        ase_gaussian.atoms = self.conformer.ase_molecule
        ase_gaussian.directory = new_scratch
        ase_gaussian.parameters["time"] = time
        ase_gaussian.parameters["partition"] = ""
        del ase_gaussian.parameters['force']

        return ase_gaussian

    def verify_output_file(self, path):
        """
        A method to verify output files and make sure that they successfully converged, if not, re-running them

        Returns a tuple where the first entry indicates if the file is complete, the second indicates if it was successful
        """

        if not os.path.exists(path):
            logging.info("Not a valid path, cannot be verified...")
            return (False, False)

        f = open(path, "r")
        file_lines = f.readlines()[-5:]
        verified = (False, False)
        for file_line in file_lines:
            if " Normal termination" in file_line:
                verified = (True, True)
            if " Error termination" in file_line:
                verified = (True, False)

        return verified

    def validate_irc(self):
        """
        A method to verify an IRC calc
        """

        logging.info("Validating IRC file...")
        irc_path = os.path.join(
            self.directory,
            "ts",
            self.conformer.reaction_label,
            "irc",
            self.conformer.reaction_label + "_irc_" + self.conformer.direction + "_" + str(self.conformer.index) + ".log"
        )

        complete, converged = self.verify_output_file(irc_path)
        if not complete:
            logging.info(
                "It seems that the IRC claculation did not complete")
            return False
        if not converged:
            logging.info("The IRC calculation did not converge...")
            return False

        pth1 = list()
        steps = list()
        with open(irc_path) as output_file:
            for line in output_file:
                line = line.strip()

                if line.startswith('Point Number:'):
                    if int(line.split()[2]) > 0:
                        if int(line.split()[-1]) == 1:
                            pt_num = int(line.split()[2])
                            pth1.append(pt_num)
                        else:
                            pass
                elif line.startswith('# OF STEPS ='):
                    num_step = int(line.split()[-1])
                    steps.append(num_step)
        # This indexes the coordinate to be used from the parsing
        if steps == []:
            logging.error('No steps taken in the IRC calculation!')
            return False
        else:
            pth1End = sum(steps[:pth1[-1]])
            # Compare the reactants and products
            irc_parse = cclib.io.ccread(irc_path)


            atomcoords = irc_parse.atomcoords
            atomnos = irc_parse.atomnos

            mol1 = rmgpy.molecule.Molecule()
            mol1.from_xyz(atomnos, atomcoords[pth1End])
            mol2 = rmgpy.molecule.Molecule()
            mol2.from_xyz(atomnos, atomcoords[-1])

            test_reaction = rmgpy.reaction.Reaction(
                reactants=mol1.split(),
                products=mol2.split(),
            )

            r, p = self.conformer.reaction_label.split("_")

            reactants = []
            products = []

            for react in r.split("+"):
                react = rmgpy.molecule.Molecule(smiles=react)
                reactants.append(react)

            for prod in p.split("+"):
                prod = rmgpy.molecule.Molecule(smiles=prod)
                products.append(prod)

            possible_reactants = []
            possible_products = []
            for reactant in reactants:
                possible_reactants.append(
                    reactant.generate_resonance_structures())

            for product in products:
                possible_products.append(
                    product.generate_resonance_structures())

            possible_reactants = list(itertools.product(*possible_reactants))
            possible_products = list(itertools.product(*possible_products))

            for possible_reactant in possible_reactants:
                reactant_list = []
                for react in possible_reactant:
                    reactant_list.append(react.to_single_bonds())

                for possible_product in possible_products:
                    product_list = []
                    for prod in possible_product:
                        product_list.append(prod.to_single_bonds())

                    target_reaction = rmgpy.reaction.Reaction(
                        reactants=list(reactant_list),
                        products=list(product_list)
                    )

                    if target_reaction.is_isomorphic(test_reaction):
                        logging.info("IRC calculation was successful!")
                        return True
            logging.info(f"IRC calculation failed for {irc_path} :(")
            return False
