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
from autotst.species import Species, Conformer
from autotst.reaction import Reaction, TS
import logging

class Orca():

    def __init__(self,conformer=None):
        
        self.command = 'orca'
        if conformer:
            assert isinstance(conformer,Conformer),'conformer must be an autotst conformer object'
            self.conformer = conformer
            if isinstance(conformer, TS):
                self.label = self.conformer.reaction_label
            else:
                self.label = self.conformer.smiles
        
    def __repr__(self):
        return '<Orca Calculator>'

    def write_fod_input(self,path):
        """
        Generates input files to run finite temperaure DFT to determine the Fractional Occupation number weighted Density (FOD number).
        Uses the default functional, basis set, and SmearTemp (TPSS, def2-TZVP, 5000 K) in Orca.
        See resource for more information:
        Bauer, C. A., Hansen, A., & Grimme, S. (2017). The Fractional Occupation Number Weighted Density as a Versatile Analysis Tool for Molecules with a Complicated Electronic Structure. 
        Chemistry - A European Journal, 23(25), 6150–6164. https://doi.org/10.1002/chem.201604682
        """
        label = self.label
        charge = self.conformer.rmg_molecule.getNetCharge()
        mult = self.conformer.rmg_molecule.multiplicity
        coords = self.conformer.get_xyz_block()

        if not os.path.exists(path):
            os.makedirs(path)
        outfile = os.path.join(path,label+'_fod.inp')

        if '(' in label or '#' in label:
            base = label.replace('(', '{').replace(')', '}').replace('#','=-')
        else:
            base = label

        with open(outfile, 'w') as f:
            f.write('# FOD anaylsis for {} \n'.format(label))
            f.write('! FOD \n')
            f.write('\n')
            f.write('%pal nprocs 4 end \n')
            f.write('%scf\n  MaxIter  600\nend\n')
            f.write('%base "{}_fod" \n'.format(base))
            f.write('*xyz {} {}\n'.format(charge, mult))
            f.write(coords)
            f.write('*\n')
    
    
    def check_NormalTermination(self,path):
        """
        checks if an Orca job terminated normally.
        Returns True is normal termination and False if something went wrong.
        """
        assert os.path.exists(path), 'It seems {} is not a valid path'.format(path)

        lines = open(path,'r').readlines()[-5:]
        complete = False
        for line in lines:
            if "ORCA TERMINATED NORMALLY" in line:
                complete = True
                break
        return complete
        
    def read_fod_log(self,path):
        """
        Reads an FOD log to get the FOD number.
        Returns FOD number if log terminated normally and FOD number can be found.
        """
        assert os.path.exists(path),'It seems {} is not a valid path'.format(path)

        if self.check_NormalTermination(path):
            N_FOD = None
            for line in open(path,'r').readlines():
                if 'N_FOD =' in line:
                    N_FOD = float(line.split(" ")[-1])
                    break
            if N_FOD:
                logging.info("the FOD number is {}".format(N_FOD))
                return N_FOD
            else:
                logging.info("It appears that the orca terminated normally for {}, but we couldn't find the FOD number".format(path))
        else:
            logging.info('It appears the orca FOD job for {} did not terminate normally'.format(path))
            

        


            



    
