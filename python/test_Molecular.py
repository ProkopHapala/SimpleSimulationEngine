#!/usr/bin/python

import numpy as np
import pyMolecular as mol

mol.lib.initWindow()




# --- spheres (atoms)
poss   = (np.random.random((100,3)) - 0.5 )*10.0
colors = np.array( [[0.9,0.1,0.1] for xyz in poss] )
radius = np.ones(len(poss))*0.5
print poss.shape, colors.shape, radius.shape
mol.spheres( poss, colors, radius )



# --- polyline (curve)
ts   = np.linspace(0,2*np.pi,100)
poss = np.transpose( np.stack([ np.sin(ts), np.cos(ts), np.sin(ts*3) ]), (1,0) ).copy()
#print poss 
mol.polyline( poss )


verts =np.array( [ [-1.0,0.0,0.0], [+1.0,0.0,0.0], [0.0,-1.0,0.0], [0.0,+1.0,0.0], [0.0,0.0,-1.0], [0.0,0.0,+1.0] ])
edges =np.array( [ [0,2], [0,3], [0,4], [0,5],   [1,2], [1,3], [1,4], [1,5],   [2,4], [2,5], [3,4], [3,5]  ], dtype=np.int32)
faces =np.array( [ [0,4,2], [0,2,5], [0,3,4], [0,5,3],   [1,2,4], [1,5,2], [1,4,3], [1,3,5] ],                dtype=np.int32)

mol.lines( edges, verts, icolor = int("0xFF10FF10", 0) )

mol.triangles( faces, verts )

mol.lib.loop(1000000)
