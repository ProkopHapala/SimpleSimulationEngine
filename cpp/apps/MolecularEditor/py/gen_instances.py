#!/usr/bin/python

import numpy as np
import sys

# === setup

fligands   = sys.argv[1]

molType = 2
L1      = 5.0
L2      = 15.0
pos0    = np.array([0.0,0.0,0.0])

prob = 0.75

# === Main

bas=np.genfromtxt( fligands, skip_header=1)

probs=np.random.rand(len(bas))
mask = probs < prob
#print( mask )


#bas=bas[:5]
#print( bas )

fout = open("instances.ini",'w')
fout.write("%i\n"  %(len(bas)+1) )
fout.write("1    %s    1.00000 0.00000 0.00000    0.00000 1.00000 0.00000  1 1\n" %(' '.join(map(str,pos0))) )
for i,l in enumerate( bas ):
	pos = l[1:]
	a   = pos - pos0
	a   = a/np.linalg.norm(a)
	b   = np.random.rand(3)
	b  -= a * np.dot(a,b)
	b   = b/np.linalg.norm(b)
	if( mask[i] ):
		pos += a*L1
	else:
		pos += a*L2
	#print( l[0], pos, a, b, np.dot(a,b) ) 
	fout.write("%i    %s    %s    %s    %i %i \n"  %( molType,' '.join(map(str,pos)),' '.join(map(str,a)),' '.join(map(str,b)), 0, 0  ) )

fout.close()

