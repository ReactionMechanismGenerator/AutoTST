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

import os, logging
from ..species import Species, Conformer
from ..reaction import Reaction, TS

class Orca():
    """
    A class for writing and reading Orca jobs.
    """

    def __init__(self, conformer=None, directory='.', partition='', time='24:00:00', nprocs=16, mem="120GB"):

        self.command = 'orca'
        self.directory = self.scratch = directory
        if conformer:
            assert isinstance(
                conformer, Conformer), 'conformer must be an autotst conformer object'
            self.conformer = conformer
            self.load_conformer_attributes()
            self.label = self.conformer.smiles
        else:
            self.label = None
            self.conformer = None

        self.nprocs = int(nprocs)
        self.mem = str(mem).upper()
        self.time = str(time)
        self.partition = partition

        if self.conformer is not None:
            num_heavy_atoms = self.conformer.rmg_molecule.get_num_atoms(
            ) - self.conformer.rmg_molecule.get_num_atoms('H')
            if num_heavy_atoms == 0:
                self.nprocs = 1
                self.mem = "10GB"
                self.time = "1:00:00"
            elif num_heavy_atoms <= 2:
                self.nprocs = 4
                self.mem = "30GB"
                self.time = "1:00:00"
            elif num_heavy_atoms <= 4:
                self.nprocs = 12
                self.mem = "60GB"
                self.time = "24:00:00"
            elif num_heavy_atoms <= 8:
                self.nprocs = 24
                self.mem = "120GB"
                self.time = "24:00:00"
            else:
                self.nprocs = 24
                self.mem = "216GB"
                self.time = "24:00:00"

        self.mem_per_proc = self.get_mem_per_proc()

    def get_mem_per_proc(self, mem=None, nprocs=None):
        """
        A method to calculate the memory per processor (in MB) for an Orca calculation.
        Returns mem_per_proc (int)

        :param mem (str or int): Total memory available in GB or MB (Assumes GB if unitless).  
        If not specified, self.mem will be used.
        :param nprocs (str or int): Number of processors.
        If not specified, self.nprocs will be used.
        """

        # Use mem is specified, else use self.mem
        if mem is None:
            assert self.mem is not None
            mem = str(self.mem).upper()
        else:
            mem = str(mem).upper()

        # Use nprocs is specified, else use self.nprocs
        if nprocs is None:
            assert self.nprocs is not None
            nprocs = int(self.nprocs)
        else:
            nprocs = int(nprocs)

        # Convert mem to MB
        if 'GB' in mem:
            mem_mb = float(mem.strip('GB')) * 1000
        elif 'MB' in mem:
            mem_mb = float(mem.strip('MB'))
        else:  # assume GB
            mem_mb = float(mem) * 1000

        # Calculate mem_per_proc
        mem_per_proc = int(mem_mb/nprocs)

        # Return mem_per_proc
        return mem_per_proc

    def __repr__(self):
        return f'<Orca Calculator {self.label}>'

    def load_conformer_attributes(self):
        """
        A method that loads attributes from a conformer attatched to an Orca calculator instance. 
        Orca calculator instance must be initialized with or provided an AutoTST conformer before calling this method.
        The method tries to get smiles, charge, multiplicity, and coordinates from the conformer.
        If found, it creates an attribute for that property.
        """

        # Assert AutoTST conformer is attached to Orca Calculator
        assert self.conformer is not None, 'Must provide an AutoTST conformer object'
        assert isinstance(
            self.conformer, Conformer), 'conformer must be an autotst conformer object'

        # Assign smiles or reaction label as label attribute
        if isinstance(self.conformer, TS):
            self.label = self.conformer.reaction_label
        else:
            self.label = self.conformer.smiles

        self.charge = self.conformer.rmg_molecule.get_net_charge()
        self.mult = self.conformer.rmg_molecule.multiplicity

        try:
            self.coords = self.conformer.get_xyz_block()
        except:
            logging.warning(
                'could not get coordinates of conformer...setting coords to None')
            self.coords = None

    def write_fod_input(self, directory=None):
        """
        Generates input files to run finite temperaure DFT to determine the Fractional Occupation number weighted Density (FOD number).
        Uses the default functional, basis set, and SmearTemp (TPSS, def2-TZVP, 5000 K) in Orca.
        See resource for more information:
        Bauer, C. A., Hansen, A., & Grimme, S. (2017). The Fractional Occupation Number Weighted Density as a Versatile Analysis Tool for Molecules with a Complicated Electronic Structure. 
        Chemistry - A European Journal, 23(25), 6150–6164. https://doi.org/10.1002/chem.201604682
        """

        # Make sure we have required properties of conformer to run the job
        assert None not in [self.mult, self.charge, self.coords]

        # If directory is not specified, use the instance directory
        if directory is None:
            directory = self.directory
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # Path for FOD input file
        outfile = os.path.join(directory, self.label+'_fod.inp')

        # Write FOD input
        with open(outfile, 'w+') as f:
            f.write(f'# FOD anaylsis for {self.label} \n')
            f.write('! FOD \n')
            f.write('\n')
            f.write(f'%pal nprocs {self.nprocs} end \n')
            f.write('%scf\n  MaxIter  600\nend\n')
            f.write('%base "fod" \n')
            f.write(f'*xyz {self.charge} {self.mult}\n')
            f.write(self.coords)
            f.write('*\n')

        return self.label+'_fod.inp'

    def write_sp_input(self, directory=None, nprocs=None, mem=None, method='ccsd(t)-f12', basis='cc-pvdz-f12',
                       scf_convergence='verytightscf', max_iter='600'):
        """
        A method to write single point energy calculation input files for ORCA.
        
        :param directory (str): directory to write input file (instance directory will be used if not specified)
        :param nprocs (str or int): number of processors to run the calculation on.
        (instance nprocs will be used if not specified)
        :param mem (str or int): amount of memory requested (in GB or MB) to run the job
        (instance mem will be used if not specified)
        :param method (str): calculation method to run.  Supported methods are (hf,ccsd(t),ccsd(t)-f12). Default is ccsd(t)-f12
        :param basis (str): basis set for calculation. Default is cc-pvdz-f12.
        :param scf_convergence : convergence option for scf. supported options are (normalscf,loosescf,sloppyscf,strongscf,tightscf,verytightscf,extremescf). Default is 'verytightscf'
        :param max_iter : maximum number of scf iterations to reach convergence criteria. (default is 600)

        Returns label for file_name
        """

        # Make sure we have required properties of conformer to run the job
        assert None not in [self.mult, self.charge, self.coords]

        self.mult = self.conformer.rmg_molecule.multiplicity
        self.charge = self.conformer.rmg_molecule.get_net_charge()
        self.coords = self.conformer.get_xyz_block()

        # Get mem_per_proc needed to write input
        if (mem or nprocs) is not None:
            if (mem and nprocs) is not None:
                mem_per_proc = self.get_mem_per_proc(mem=mem, nprocs=nprocs)
            elif nprocs is not None:
                mem = self.mem
                mem_per_proc = self.get_mem_per_proc(nprocs=nprocs)
            else:
                nprocs = self.nprocs
                mem_per_proc = self.get_mem_per_proc(mem=mem)
        else:
            mem_per_proc = self.mem_per_proc
            mem = self.mem
            nprocs = self.nprocs

        # if directory is not specified,
        # Use Orca instance directory
        if directory is None:
            directory = self.directory
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Convert parameters to lowercase strings
        method = method.lower()
        basis = basis.lower()
        scf_convergence = scf_convergence.lower()
        max_iter = max_iter.lower()

        # Assert scf_convergence is supported
        scf_convergence_options = 'normalscf loosescf sloppyscf strongscf tightscf verytightscf extremescf'
        assert scf_convergence in scf_convergence_options

        # Select UHF for open-shell,
        # RHF for closed shell molecules
        if 'hf' not in method:
            if int(self.conformer.rmg_molecule.multiplicity) == 1:
                hf = 'RHF'
            else:
                hf = 'UHF'
        else:
            if int(self.conformer.rmg_molecule.multiplicity) == 1:
                method = 'RHF'
                hf = ''
            else:
                method = 'UHF'
                hf = ''

        # Auxiliary basis sets for F12 methods
        auxiliary_basis_sets_dict = {
            'cc-pvdz-f12': 'CC-PVDZ-F12-CABS CC-PVTZ/C',
            'cc-pvtz-f12': 'CC-PVTZ-F12-CABS CC-PVQZ/C',
            'cc-pvqz-f12': 'CC-PVQZ-F12-CABS CC-PVQZ/C'}

        # Get auxiliary basis sets
        if basis in auxiliary_basis_sets_dict.keys():
            aux_basis = auxiliary_basis_sets_dict.get(basis)
        else:
            aux_basis = ''

        # Create file name
        file_name = f"{self.label}_{method}_{basis}.inp"

        # Create label
        label = file_name.strip('.inp')

        # Create file path
        file_path = os.path.join(directory, file_name)

        # Create base name for scratch files
        base = method
        if '(' in base or '#' in base or '/' in base:
            base = base.replace('(', '{').replace(')', '}').replace('#', '=-').replace('/','-')

        # Write sp input
        with open(file_path, 'w') as f:
            f.write(
                f'# {method}/{basis} calculation for {self.label}\n')
            if "f12" in method:
                f.write(
                    f'! {method}/RI {basis} {aux_basis} {scf_convergence} FrozenCore UNO UCO\n')
            else:
                f.write(
                    f'! {method} {basis} {aux_basis} {scf_convergence} FrozenCore UNO UCO\n')
            f.write(f'%pal nprocs {nprocs} end \n')
            f.write(f'%maxcore {mem_per_proc}\n')
            f.write(f'%mdci\n  MaxIter  {max_iter}\nend\n')
            f.write(f'%base "{base}"\n')
            f.write(f'*xyz {self.charge} {self.mult}\n')
            f.write(self.coords)
            f.write('*\n')

        return label

    def read_energy(self, path):
        """
        Returns the final single point energy (in Hartree) from an Orca calculation
        Checks if the job terminated normally; if not, returns False
        """

        if self.check_normal_termination(path) is False:
            logging.warning(f"It appears {path} orca calculation failed")
            return False

        lines = open(path, 'r').readlines()
        for line in lines:
            if "FINAL SINGLE POINT ENERGY" in line:
                return float(line.split()[-1].split('\n')[0])
        

    def check_normal_termination(self, path):
        """
        Checks if an Orca job terminated normally.
        Returns True is normal termination and False if something went wrong.
        """
        assert os.path.exists(
            path), 'It seems {} is not a valid path'.format(path)

        lines = open(path, 'r').readlines()[-5:]
        complete = False
        for line in lines:
            if "ORCA TERMINATED NORMALLY" in line:
                complete = True
                break
        if complete is True:
            logging.info(f"{path} completed successfully")
        else:
            logging.info(f"{path} failed")
        return complete

    def read_fod_log(self, path):
        """
        Reads an FOD log to get the FOD number.
        Returns FOD number if log terminated normally and FOD number can be found.
        """
        assert os.path.exists(path),f'It seems {path} is not a valid path'

        if self.check_normal_termination(path):
            N_FOD = None
            for line in open(path, 'r').readlines():
                if 'N_FOD =' in line:
                    N_FOD = float(line.split(" ")[-1])
                    break
            if N_FOD:
                logging.info(f"the FOD number is {N_FOD}")
                return N_FOD
            else:
                logging.info(
                    f"It appears that the orca terminated normally for {path}, but we couldn't find the FOD number")
        else:
            logging.info(
                f'It appears the orca FOD job for {path} did not terminate normally')
