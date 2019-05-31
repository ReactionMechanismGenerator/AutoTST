# Automated Transition State Theory Calculator | AutoTST
[![Build Status](https://travis-ci.com/ReactionMechanismGenerator/AutoTST.svg?branch=master)](https://travis-ci.com/ReactionMechanismGenerator/AutoTST)
[![codecov](https://codecov.io/gh/reactionmechanismgenerator/autotst/branch/master/graph/badge.svg)](https://codecov.io/gh/reactionmechanismgenerator/autotst)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
![Release](https://img.shields.io/badge/version-0.1.0-brightgreen.svg)

## Descritpion

AutoTST is a framework to perform automated transition state theory calculations related to reaction families common in combustion.
It is based off of the frame work of [**RMG**](rmg.mit.edu), [**RDKit**](http://www.rdkit.org/), and [**ASE**](https://wiki.fysik.dtu.dk/ase/) to automatically create three dimentional geometries of transition states, reactants, and products;
optimize them using the [**Gaussian 09**](http://gaussian.com/) quantum chemistry package; and obtain kinetic parameters using **Arkane** (a canonical calculator included in RMG).

Currently, AutoTST supports three reaction families:
- Hydrogen Abstraction
- Disproportination
- Intra Hydrogen Migration

For general templates of these reaction families, click [here](https://github.com/ReactionMechanismGenerator/RMG-database/blob/master/families/rmg_reaction_families.pdf)

However, we intend to introduce support for new reaction families and quantum packages.


## How to Install

Before installing AutoTST, download [Anaconda](anaconda.com/download/) and [Git](https://git-scm.com/downloads)

Install the latest version of AutoTST by cloning the source code via Git. Make sure to start in an appropriate local directory where you want the AutoTST folder to exist.

- `git clone https://github.com/ReactionMechanismGenerator/AutoTST.git`

Now, create the anaconda environment for AutoTST

- `cd AutoTST`

- `conda env create -f environment.yml`

Finally, RMG has a dependency on `openbabel` but this dependency causes AutoTST to crash. So, we have to uninstall it with the following line of code:

- `conda uninstall --force openbabel`


Modify environment variables. Add AutoTST to the `PYTHONPATH` to ensure that you can access modules from any folder. Modify your ~/.bashrc file by adding the following line:

- `export PYTHONPATH=$PYTHONPATH:your_folder/AutoTST`

To be able to run AutoTST in any conda environment, you can set your path to the following by modifing your ~/.bashrc: 

- `export PATH=~/anaconda/envs/tst_env/bin:$PATH`

Finally either close and reopen your terminal to refresh your environment variables, or type the following command 
- `source ~/.bashrc`

## How to Give Feedback

Please post any issues you may have to the [issues page](https://github.com/ReactionMechanismGenerator/AutoTST/issues/)
or drop in to the [chat room](https://gitter.im/ReactionMechanismGenerator/AutoTST) 

## Credits

- [Professor Richard H. West's research group](http://www.northeastern.edu/comocheng/) at 
[Northeastern University](http://www.northeastern.edu/). 
- Dr. Pierre L. Bhoorasingh
