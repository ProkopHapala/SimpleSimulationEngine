#!/usr/bin/python

import numpy as np
import pyMolecular as mol

mol.initWorld("/home/prokop/git/SimpleSimulationEngine/cpp/apps/MolecularEditor/inputs/")
mol.relax( 1500 )
mol.exportAtoms( "final.xyz" )

