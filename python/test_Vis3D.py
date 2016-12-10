#!/usr/bin/python

import numpy as np
import pyVis3D as vis
from time import sleep

vis.lib.initWindow()
vis.lib.asyncLoop(1000000000)  # asynchronous window - does not wait for exit

# --- spheres (atoms)
poss_sph = (np.random.random((100,3)) - 0.5 )*10.0
colors   = np.array( [[0.9,0.1,0.1] for xyz in poss_sph] )
radius   = np.ones(len(poss_sph))*0.5
print poss_sph.shape, colors.shape, radius.shape
glint_spheres = vis.spheres( poss_sph, colors, radius )


# --- polyline (curve)
ts   = np.linspace(0,2*np.pi,100)
poss = np.transpose( np.stack([ np.sin(ts), np.cos(ts), np.sin(ts*3) ]), (1,0) ).copy()
#print poss 
glint_curve = vis.polyline( poss )

verts =np.array( [ [-1.0,0.0,0.0], [+1.0,0.0,0.0], [0.0,-1.0,0.0], [0.0,+1.0,0.0], [0.0,0.0,-1.0], [0.0,0.0,+1.0] ])
edges =np.array( [ [0,2], [0,3], [0,4], [0,5],   [1,2], [1,3], [1,4], [1,5],   [2,4], [2,5], [3,4], [3,5]  ], dtype=np.int32)
faces =np.array( [ [0,4,2], [0,2,5], [0,3,4], [0,5,3],   [1,2,4], [1,5,2], [1,4,3], [1,3,5] ],                dtype=np.int32)

glint_lines = vis.lines( edges, verts, icolor = int("0xFF10FF10", 0) )

glint_trinagles = vis.triangles( faces, verts )

#vis.lib.loop(1000000)  # synchronous window - wait for exit

for i in range(20):
    sleep(0.5)
    print( glint_curve, glint_spheres );
    # --- change curve shape
    vis.lib.erase( glint_curve )
    poss        = np.transpose( np.stack([ np.sin(ts), np.cos(ts), np.sin(ts*i) ]), (1,0) ).copy()
    glint_curve = vis.polyline( poss )
    # --- change sphere colors
    vis.lib.erase( glint_spheres );                                    print "DEBUG 1 "
    colors = np.array( [[0.05*i,0.5,1-0.05*i] for xyz in poss_sph] );  print "DEBUG 2 "
    glint_spheres = vis.spheres( poss_sph, colors, radius );           print "DEBUG 3 "

'''
vis.lib.initWindow()
vis.lib.asyncLoop(10000)
for i in range(10):
    vis.lib.setGobVar( i*100.0001 );
    sleep(0.05)
'''




