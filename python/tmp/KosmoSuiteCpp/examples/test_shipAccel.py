#!/usr/bin/python

from pylab import *
import numpy as np
import os 


print "========= Load Library ============="

def makeclean():
	import os
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".so") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".o") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".pyc") ]

os.chdir('../');       print " >> WORKDIR: ", os.getcwd()
makeclean( )
sys.path.insert(0, "./")
import KosmoSuiteCpp as KS
os.chdir('examples');  print " >> WORKDIR: ", os.getcwd()

print "========= Setup ============="

trj_tin  = linspace( 0, 200000, 1000 );   nstep = len(trj_tin) 
trj_tout = zeros( nstep ) 
trj_v    = zeros( ( nstep, 3 ) );
trj_r    = zeros( ( nstep, 3 ) );

r0 = array( [ 6771.0e+3, 0.0,    0.0 ] ) 
v0 = array( [    0.0,    7.9e+3, 0.0 ] ) 
CM = 5.97219e+24 * 6.67384e-11;


accels = [  1e-2, 2e-2, 4e-2, 4e-2 ]

figure()
for accel in accels:
	KS.shipAccel_setup( -CM, accel, 1e-3, 1e-5, r0, v0 )
	KS.shipAccel_run  ( 100.0, 10.0, 1000.0, trj_tin, trj_tout, trj_r, trj_v )
	print trj_r
	#plot( trj_r[:,0].copy(), trj_r[:,1].copy(), '.-', label='%e' %accel )
	plot( trj_r[:,0].copy(), trj_r[:,1].copy(), '-', label='a=%2.2e m/s2' %accel )

legend()
axis('equal')
grid()

show()

print " all DONE ! "


 
