#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from cclib.io import ccread

__file_path = os.path.dirname(os.path.abspath(__file__))

settings = {'autotst_path': __file_path,
            'tst_database_path': os.path.join(os.path.dirname(__file_path), 'database'),
            }

