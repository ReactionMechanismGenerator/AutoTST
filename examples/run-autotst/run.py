import time
import os
import pickle
import subprocess
import multiprocessing
from multiprocessing import Process
import logging
import subprocess
import argparse
FORMAT = "%(filename)s:%(lineno)d %(funcName)s %(levelname)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.INFO)
from autotst.job.job import *
from autotst.calculator.gaussian import Gaussian
from ase.calculators.dftb import Dftb
from multiprocessing import Manager, Process
from rmgpy.molecule import Molecule as RMGMolecule 
from rmgpy.species import Species as RMGSpecies
from rmgpy.reaction import Reaction as RMGReaction

parser = argparse.ArgumentParser(description='Command Line options for running autotst')

parser.add_argument('--calculation-type', default='thermo', help='thermo or kinetics calculation')
parser.add_argument('--scracth-dir',default=os.getenv('AUTOTST_SCR_DIR'), help='scracth directory to store your files')
parser.add_argument('--reaction-family', default='H_Abstraction', help='Reaction family for autotst')
parser.add_argument('--partition', default='short', help="partition to run your calculations")
parser.add_argument('--autotst-label', default='[O]O', help='Smiles string for thermo or reaction string for kinetics' )

def main():
    args = parser.parse_args()


    if os.getenv('SLURM_ARRAY_TASK_ID'):
        i = int(os.getenv('SLURM_ARRAY_TASK_ID'))
    else:
        #raise Exception("Specif y a TS number!")
        logging.warning("Number not specified as script argument or via environment variable, so using default")
        i = 5
    #i = i + 999 #### ADDED THIS LINE TO GET ARRAY OVER 1000 FOR SLURM
    logging.info("RUNNING WITH JOB NUMBER i = {}".format(i))
    logging.info("Running on a node with {} processors".format(multiprocessing.cpu_count()))
    
    if args.calculation_type == 'thermo':
        autotst_species = Species(smiles=[args.autotst_label])
        job = Job(
        calculator=Gaussian(
            settings={
                "method": "m062x",
                "basis": "cc-pVTZ",
                "mem": "5GB",
                "nprocshared": 20,
            },
            directory = args.scracth_dir 
        ),
        username=os.getenv("USER"),
        conformer_calculator=Dftb(directory=args.scracth_dir),
        partition=args.partition,
        scratch=args.scracth_dir
        )
        logging.info("Submitting low energy conformers for thermo calculation")
        result = job.calculate_species(species=autotst_species)
    elif args.calculation_type == 'kinetics':
        #autotst_reaction= Reaction(label=args.autotst_label,reaction_family=args.reaction_family)
        reaction = Reaction(label=args.autotst_label, reaction_family=args.reaction_family)
        job = Job(
        reaction=reaction,
        calculator=Gaussian(
            settings={
                "method": "m062x",
                "basis": "cc-pVTZ",
                "mem": "5GB",
                "nprocshared": 20,
            },
            directory = args.scracth_dir
        ),
        username=os.getenv("USER"),
        conformer_calculator = Dftb(directory=args.scracth_dir),
        partition=args.partition,
        scratch=args.scracth_dir
        )
        logging.info("Submitting transition states for kinetics calculation")
        result = job.calculate_reaction()
        

    if result:
        print("AutoTST successfully arrived at an optimized geometry for {}".format(autotst_species))
    else:
        print("AutoTST could not arrive at an optimized geometry for {} :(".format(autotst_species))

    print("All geometries completed")
    print("Done!")

if __name__ == '__main__':
    main()