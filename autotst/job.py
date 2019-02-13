import numpy as np
import logging
import subprocess
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)

import rdkit, rdkit.Chem.rdDistGeom, rdkit.DistanceGeometry

from rdkit import Chem

from rdkit.Chem import AllChem
from rdkit.Chem.Pharm3D import EmbedLib

import ase

import rmgpy
from rmgpy.molecule import Molecule as RMGMolecule
from rmgpy.species import Species as RMGSpecies
from rmgpy.reaction import Reaction as RMGReaction, ReactionError
from rmgpy.kinetics import PDepArrhenius, PDepKineticsModel
from rmgpy.data.rmg import RMGDatabase

import os
from autotst.reaction import Reaction, TS
from autotst.molecule import Molecule
from autotst.geometry import Bond, Angle, Torsion, CisTrans, ChiralCenter
from ase.io.gaussian import read_gaussian_out

# To perform TS search
from autotst.calculators.calculator import Calculator
from ase.calculators.gaussian import Gaussian as ASEGaussian
from autotst.calculators.gaussian import Gaussian
from autotst.calculators.vibrational_analysis import VibrationalAnalysis, percent_change
from autotst.calculators.statmech import StatMech

 

class Job():
    """
    A class to deal with the input and output of calculations
    """

    def __init__(self, reaction=None, dft_calculator=None, conformer_calculator=None):

        self.reaction = reaction
        if isinstance(reaction, Reaction):
            self.label = reaction.label

        self.dft_calculator = dft_calculator
        self.conformer_calculator = conformer_calculator

    def __repr__(self):
        return "< Job '{}'>".format(self.calculator, self.label)


    def submit(self, conformer=None, ase_calculator=None, partition="general"):

        if not ase_calculator:
            logging.info("Cannot submit job... No calculator provided.")
            return None

        self.write_input(conformer, ase_calculator)

        command = ase_calculator.command
        label = ase_calculator.label

        submit_string = "srun -p {0} -N 1 -n 16 --mem=128GB {1} < {2}.com > {2}.log".format(
            partition,
            command,
            label
        )

        errorcode = subprocesses.call(submit_string, shel=True, cwd=ase_calculator.scratch)

        if errorcode:
            logging.error("SOMETHING WENT WRONG")
            return None

        else:
            return "Yay"


    def write_input(self, conformer=None, ase_calculator=None):
        """
        A method to write the IO file needed for a quantum calculation
        """
        if not conformer:
            assert ase_calculator.atoms, "Calculator needs to have atoms object attached"

        else:
            ase_calculator.atoms = conformer.ase_molecule

        ase_calculator.write_input()


    def run_rotors(self, calculator=None):

        return None
    
    def run_conformer_analysis(self, species=None, calculator=None):
        """
        A method to run conformer analysis at a low level of theory
        """
        from hotbit import Hotbit() # only using hotbit for now
        conformers = species.generate_conformers(method="systematic", calculator=Hotbit())

        complete = []
        for smiles, confs in conformers.iteritems():
            for conf in confs:
                calc = calculator.get_species_calc(conf)
                self.submit(conf, calc)
                complete.append(False)

        from cclib.io import ccread
        while not np.array(complete).all():
            i = 0
            for smiles, confs in conformers.iteritems():
                for conf in confs:
                    calc = calculator.get_species_calc(conf)
                    p = ccread(os.path.join(calc.scratch, calc.label +".log"))
                    if p.optdone == True:
                        complete[i] = True
                    i += 1


        # assuming that all quantum calcs are done in gaussian

        return None

    def run_species_optimization(self, conformer=None, calculator=None):
        """
        A method for optimizing geometries for species
        """

        assert (isinstance(conformer, Conformer) and (not isinstance(conformer, TS))), "Please provide a Conformer object."
        assert calculator
            
            calc = calculator.get_species_calc(conformer)
            conformer, result = calculator.calculate(conformer, calc)
            self.fix_io_file(calc)
            
            if result:
                logging.info("TS validated, now running hindered rotor calculations")
                ### Add hindered rotor work here
                logging.info("jk, this feature hasn't been added just yet")
                
            if result:
                logging.info("Conformer species successfully optimized")
                return result
            
            else:
                logging.info("Could not optimize species geometry")
                return result

    def run_ts_search(self, conformer=None, calculator=None, vibrational_analysis=True, hindered_rotors=True):
        """
        A method to run a transition state search
        """
        
        assert isinstance(conformer, TS), "`conformer` provided not a TS type..."
        
        # Performing the TS optimizations
        logging.info("Conformer provided is a TS object")
        
        shell = calculator.get_shell_calc(conformer)
        logging.info("Running optimization of reaction shell")
        conformer, result = calculator.calculate(conformer, shell)
        calculator.fix_io_file(shell)
        if not result:
            logging.info("FAILED SHELL CALCULATION")
            return result
        
        center = calculator.get_center_calc(conformer)
        logging.info("Running optization of reaction center")
        conformer, result = calculator.calculate(conformer, center)
        calculator.fix_io_file(center)
        if not result:
            logging.info("FAILED CENTER CALCULATION")
            return result
        
        overall = calculator.get_overall_calc(conformer)
        logging.info("Running overall optimization of TS")
        conformer, result = calculator.calculate(conformer, overall)
        calculator.fix_io_file(overall)
        if not result:
            logging.info("FAILED OVERALL CALCULATION")
            
        if not vibrational_analysis:
            logging.info("Running without vibrational analysis. \nRunning IRC instead")
            irc = calculator.get_irc_calc(conformer)
            calculator.run_irc(conformer, irc)
            result = calculator.validate_irc(irc)
            calculator.fix_io_file(irc)

        else:
            from autotst.calculators.vibrational_analysis import VibrationalAnalysis
            vib = VibrationalAnalysis(ts=conformer, scratch=calculator.scratch)
            result = vib.validate_ts()
            
            if not result:
                logging.info("Vibrational Analysis not conclusive...\n Running IRC instead")
                irc = calculator.get_irc_calc(conformer)
                calculator.run_irc(conformer, irc)
                result = calculator.validate_irc(irc)
                calculator.fix_io_file(irc)
                
        if (result and hindered_rotors):
            logging.info("TS validated, now running hindered rotor calculations")
            ### Add hindered rotor work here
            logging.info("jk, this feature hasn't been added just yet")
                
        if result:
            logging.info("Arrived at a TS!")
            return result

        else:
            logging.info("Could not arrive at a TS!")
            return result


    def run(self, 
        conformer=None, 
        vibrational_analysis=True, 
        hindered_rotors=True,):
        """
        A method to perform all the necessary calculations required for a particular conformer
        """
        
        if not conformer:
            conformer = self.conformer
        
        assert isinstance(conformer, (Conformer, TS)), "`conformer` provided not a Conformer type..."
        
        if isinstance(conformer, TS):
            # Performing the TS optimizations
            logging.info("Conformer provided is a TS object")
            
            shell = self.get_shell_calc(conformer)
            logging.info("Running optimization of reaction shell")
            conformer, result = self.calculate(conformer, shell)
            self.fix_io_file(shell)
            if not result:
                logging.info("FAILED SHELL CALCULATION")
                return result
            
            center = self.get_center_calc(conformer)
            logging.info("Running optization of reaction center")
            conformer, result = self.calculate(conformer, center)
            self.fix_io_file(center)
            if not result:
                logging.info("FAILED CENTER CALCULATION")
                return result
            
            overall = self.get_overall_calc(conformer)
            logging.info("Running overall optimization of TS")
            conformer, result = self.calculate(conformer, overall)
            self.fix_io_file(overall)
            if not result:
                logging.info("FAILED OVERALL CALCULATION")
                
            if not vibrational_analysis:
                logging.info("Running without vibrational analysis. \nRunning IRC instead")
                irc = self.get_irc_calc(conformer)
                self.run_irc(conformer, irc)
                result = self.validate_irc(irc)
                self.fix_io_file(irc)

            else:
                from autotst.calculators.vibrational_analysis import VibrationalAnalysis
                vib = VibrationalAnalysis(ts=conformer, scratch=self.scratch)
                result = vib.validate_ts()
                
                if not result:
                    logging.info("Vibrational Analysis not conclusive...\n Running IRC instead")
                    irc = self.get_irc_calc(conformer)
                    self.run_irc(conformer, irc)
                    result = self.validate_irc(irc)
                    self.fix_io_file(irc)
                    
            if result:
                logging.info("TS validated, now running hindered rotor calculations")
                ### Add hindered rotor work here
                logging.info("jk, this feature hasn't been added just yet")
                    
            if result:
                logging.info("Arrived at a TS!")
                return result

            else:
                logging.info("Could not arrive at a TS!")
                return result
            
        elif isinstance(conformer, Conformer):
            logging.info("Conformer provided is NOT a TS object")
            
            calc = self.get_species_calc(conformer)
            conformer, result = self.calculate(conformer, calc)
            self.fix_io_file(calc)
            
            if result:
                logging.info("TS validated, now running hindered rotor calculations")
                ### Add hindered rotor work here
                logging.info("jk, this feature hasn't been added just yet")
                
            if result:
                logging.info("Conformer species successfully optimized")
                return result
            
            else:
                logging.info("Could not optimize species geometry")
                return result

    

    

    