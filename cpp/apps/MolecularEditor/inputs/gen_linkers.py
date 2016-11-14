#!/usr/bin/python

import numpy as np

#fin = open("H56.bas",r)

bas=np.genfromtxt("H56.bas", skip_header=1)

#bas=bas[:5]
#print( bas )


k0=30.0
l0=0.0
mol0=0
pos0 = np.array([0.0,0.0,0.0])

fout = open("linkers.ini",'w')
fout.write("%i\n"  %(len(bas)) )

for i,l in enumerate( bas ):
	pos = l[1:]
	fout.write("%i %i    %s  %s   %f  %f\n"  %( mol0,i+1, ' '.join(map(str,pos)),' '.join(map(str,pos0)),   k0, l0 ) )

fout.close()

