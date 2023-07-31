"""

This module contains common functions that are used throughout AutoTST
Therefore, it should not import any other AutoTST modules to avoid circular imports

"""

import os

import rmgpy


# Absolute path to the AutoTST directory
AUTOTST_PATH = os.path.dirname(os.path.abspath(__file__))
# Abolsute path to RMG-Py directory
RMG_PY_PATH = os.path.abspath(os.path.dirname(os.path.dirname(rmgpy.__file__)))
# Absolute path to the RMG-database directory
RMG_DATABASE_PATH = os.path.abspath(os.path.dirname(os.path.dirname(rmgpy.settings['database.directory'])))

VERSION = "1.0.0"
