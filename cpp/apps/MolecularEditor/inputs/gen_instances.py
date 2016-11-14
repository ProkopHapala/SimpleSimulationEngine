#!/usr/bin/python

import numpy as np

#fin = open("H56.bas",r)

bas=np.genfromtxt("H56.bas", skip_header=1)

#bas=bas[:5]
#print( bas )

molType=2

L0 = 5.0

fout = open("instances.ini",'w')
fout.write("%i\n"  %(len(bas)+1) )
fout.write("1    0.00000  0.00000  0.00000    1.00000 0.00000 0.00000    0.00000 1.00000 0.00000  1 1\n" )
for i,l in enumerate( bas ):
	pos = l[1:]
	a   = pos.copy()
	a   = a/np.linalg.norm(a)
	b   = np.random.rand(3)
	b  -= a * np.dot(a,b)
	b   = b/np.linalg.norm(b)

	pos += a*L0
	#print( pos, a, b, np.dot(a,b) ) 
	fout.write("%i    %s    %s    %s    %i %i \n"  %( molType,' '.join(map(str,pos)),' '.join(map(str,a)),' '.join(map(str,b)), 0, 0  ) )

fout.close()

