#!/usr/bin/python

from pylab import *
import numpy as np
from pySimE.KosmoSuiteCpp import main as KS

print " =========== sample magnetic field ============="

n = 100
ps  = zeros((n,3))
dIs = zeros((n,3))

m = 100
wheres = zeros((m,3))
Bs     = zeros((m,3))

R  = 1.0 
dl = 2*pi*R/n;
I  = 1.0e+6; 
#phis = linspace( 0, 2*pi, n ) 
phis = arange( 0, 2*pi, 2*pi/n ) 
ca   = cos(phis)
sa   = sin(phis)
ps [:,0] = ca * R 
ps [:,1] = sa * R
dIs[:,0] =  sa;
dIs[:,1] = -ca;
dIs   *= dl*I; 

#wheres[:,0] = 0.0
#wheres[:,1] = 0.9
#wheres[:,2] = linspace(-5,5, m) 

wheres[:,0] = 0.0 
wheres[:,1] = linspace(-5,5, m)
wheres[:,2] = 1.0 

#print ps
#print dIs
#print wheres

KS.elmag_sample( wheres, Bs, ps, dIs )

#print Bs

plot( wheres[:,1], Bs[:,0] )
plot( wheres[:,1], Bs[:,1] )
plot( wheres[:,1], Bs[:,2] )

grid()

print " =========== particle trajectory ============="

eVtoSI = 1.60217657e-19

charge  = 1.60217656535e-19; # [C]
mass    = 1.67262178e-27;    # [kg]   me=9.10938291e-31, mp = 1.67262178e-27   ma = 1.66053892173e-27
E       = 1e+5;
ESI     = E *  eVtoSI;
vSI     = sqrt( ( 2 * ESI ) / mass ) 
print " ESI, vSI ", ESI, vSI

trj_tin  = linspace( 0, 0.2e-5, 1000 );   nstep = len(trj_tin) 
trj_tout = zeros( nstep ) 
trj_v    = zeros( ( nstep, 3 ) );
trj_r    = zeros( ( nstep, 3 ) );

#figure();  plot( ps[:,0], ps[:,1], '.-k', label='coil' )
figure();  plot( ps[:,0], ps[:,1], '-k', lw=2, label='coil' )

r0 = array( [   -3, 0.0, 1.0 ] ) 
v0 = array( [  vSI, 0.0, 0.0 ] ) 
ys = linspace( -1.0, 1.0, 10 )
for y in ys:
	r0[1] = y;
	KS.elmag_setup( charge, mass, 1e-9, 1e+0, r0, v0, ps, dIs )
	KS.elmag_run  ( 1.0e-9, 1.0e-10, 1.0e-8, trj_tin, trj_tout, trj_r, trj_v )
	#plot( trj_r[:,0], trj_r[:,1], '.-', label='%2.2e' %y )
	plot( trj_r[:,0], trj_r[:,1], '-', label='y=%2.2f m' %y )

legend(loc=2)
grid()
axis('equal')

show()

 
