# -*- coding: utf-8 -*-
"""
Created on Sat May 31 14:59:04 2014

@author: asiJa
"""

import sys
bak_modules = sys.modules.copy()

def reset_sys():
	for k in sys.modules.keys():
	    if not k in bak_modules:
	        del sys.modules[k]	






