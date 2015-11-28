#!/usr/bin/python

from pylab import *
import numpy as np
from pySimE.KosmoSuiteCpp import main as KS

# ============ Setup variables

Rs   = array([ [0.0,0.0,0.0],  [149.598261e+9,       0.0, 0.0 ] ] ).copy()   
Vs   = array([ [0.0,0.0,0.0],  [          0.0, 29.780e+3, 0.0 ] ] ).copy()
Ms   = array([ 1.9891e+30,       5.97219e+24                  ]   ).copy()
Errs = array([ [1e-1,1e-5],    [ 1e-1, 1e-5                   ] ] ).copy()
nbodies = len(Ms)

print "nbodies ", nbodies 

trj_tin  = linspace( 0, 10000000, 100 );   nstep = len(trj_tin) 
#trj_tin  = linspace( 0, 1000, 2 );   nstep = len(trj_tin) 
print trj_tin
trj_tout = zeros( nstep ) 
trj_v    = zeros( ( nstep, nbodies, 3 ) );
trj_r    = zeros( ( nstep, nbodies, 3 ) );

# ============ Simulation

print " KS.nbody_setup( Ms, Rs, Vs ) "
KS.nbody_setup( Ms, Rs, Vs, Errs )
print " KS.nbody_run  ( 1000.0, trj_tin, trj_tout, trj_r, trj_v ) "
KS.nbody_run  ( 10000.0, 1000.0, 100000.0, trj_tin, trj_tout, trj_r, trj_v )

print trj_r

print " simulation done ! "

figure()
plot( trj_tin, trj_tout, '.-' )

figure();
plot( trj_r[:,1,0], trj_r[:,1,1], '.-' )
axis('equal')


show()


print " all DONE ! "


 
