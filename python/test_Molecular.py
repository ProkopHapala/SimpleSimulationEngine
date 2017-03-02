#!/usr/bin/python

import numpy as np
import pyMolecular as mol

mol.initWorld("/home/prokop/git/SimpleSimulationEngine/cpp/apps/MolecularEditor/inputs/")

# debug - test output I/O for molecular geometry
arr = mol.getInstancePointer()
print "opt array ", arr
types, poss = mol.getAtoms()
print len(types)
print types 
print poss


mol.relax( 1500 )
mol.exportAtoms( "final.xyz" )

