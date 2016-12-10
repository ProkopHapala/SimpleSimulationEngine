#!/usr/bin/python

import numpy as np
import pyMolecular as mol

points_ref = np.array([
 [1.0,0.0,0.0],   [-1.0, 0.0, 0.0], 
 [0.0,1.0,0.0],   [ 0.0,-1.0, 0.0],
 [0.0,0.0,1.0],   [ 0.0, 0.0,-1.0]
], dtype=np.float64 )

mol.initComparator( points_ref )

points = points_ref.copy()
dist   = mol.compDistance( points );   print( "dist (identical)", dist )

drnd    = np.random.rand( points.shape[0], points.shape[1] )
points += drnd*0.01
print "========== drandom"
print points
dist   = mol.compDistance( points.copy() );   print( "dist (drandom)", dist )

np.random.shuffle(points)
print points
dist = mol.compDistance( points.copy() );   print( "dist (shuffled)", dist )
