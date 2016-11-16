#!/usr/bin/python

import numpy as np
import sys

# === setup

k0   = 50.0   # this converge with time_step 0.01 
#k0    = 5.0  # this converge with time_step 0.05
l0    = 1.94
dlmax = 0.5
imol  = 0
jatom = 0
pos0 = np.array([0.0,0.0,0.0])

fligands   = sys.argv[1]
fsubstrate = sys.argv[2]

# === Functions

def findNearest( p, ps ):
	rs = np.sum( (ps - p)**2, axis=1 )
	#print (rs)
	return np.argmin(rs)

# === Main

bas =np.genfromtxt( fligands,   skip_header=1)
bas2=np.genfromtxt( fsubstrate, skip_header=1)

#bas=bas[:5]

poss =bas2[:,1:]

fout = open("bonds.ini",'w')
fout.write("%i\n"  %(len(bas)) )

fDEBUG = open("DEBUG.bas",'w')
fDEBUG.write("%i\n"  %(len(bas2) + len(bas)) )

for jmol,l in enumerate( bas ):
	pos = l[1:]
	iatom = findNearest( pos, poss )
	fout.write("%i %i    %i  %i   %f %f %f\n"  %( imol, jmol+1,   iatom, jatom,   k0, l0, dlmax ) )

	bas2[iatom,0] = 20
	fDEBUG.write("%i %f %f %f\n"  %( 1, pos[0], pos[1], pos[2] ) ) 


for i,l in enumerate( bas2 ):
	fDEBUG.write("%i %f %f %f\n"  %( l[0], l[1], l[2], l[3] ) ) 
fDEBUG.close()

fout.close()

