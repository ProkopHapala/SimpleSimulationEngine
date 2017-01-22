#!/usr/bin/python
# use like:
# python replace.py answer.bas 1 14 2.0 0.75


import sys
import numpy as np
import atomicUtils as au
import plotUtils as pu
import matplotlib.pyplot as plt

atoms          = np.genfromtxt  ( 'answer.bas', skip_header=1 )
bonds,bondVecs = au.findAllBonds( atoms, Rcut=3.0, RvdwCut=0.8 )
neighs         = au.neighs( len(atoms), bonds )

select1 = au.findTypeNeigh( atoms, neighs, 14, neighTyps={1:(2,2)} )

select2 = au.getAllNeighsOfSelected( select1, neighs, atoms, typs={1} )
select2 = list(select2.keys())

#atoms[select1,0] = 16
#atoms[select2,0] = 9
#au.saveAtoms( atoms, "test.xyz", xyz=True )

group = np.array([
[6, 0.0, 0.0, 0.0 ],
[1, +0.7,0.7, 0.0 ],
[1, -0.7,0.7, 0.0 ]
])

group = np.array([
[15, 0.0, 0.0, 0.0 ],
[1, 0.0,1.5, 0.0 ],
])


'''
group = np.array([
[6, 0.0, 0.0, 0.0 ],
[8, 1.0, 0.0, 0.0 ],
[9, 0.0, 1.0, 0.0 ],
[7, 0.0, 0.0, 1.0 ]
])
'''

pairs  = au.findPairs_one( select2, atoms, Rcut=2.5 );   #print( pairs )

pairs = au.pairsNotShareNeigh( pairs, neighs );   #print( pairs )

cog = au.findCOG( atoms[:,1:] )

atoms_ = au.replacePairs( pairs, atoms, group, up_vec=(cog,1) );   # print( "atoms_ = ",atoms_ )
au.saveAtoms( atoms_, "atoms_.xyz", xyz=True )

# ---------- ploting

'''
rotMat = au.makeRotMat( [1.0,1.0,1.0], [0.0,1.0,0.0] )
ps     = np.dot( rotMat, np.transpose(atoms[:,1:]) )


print( ps.shape )

#plt.plot( ps[0], ps[1], 'ok' )

#pu.plotAtoms( atoms[:,0], ps[0], ps[1] )
#pu.plotBonds( [ b[0] for b in bonds], ps[0], ps[1] )
#pu.plotBonds( pairs, ps[0], ps[1] )

plt.axis('equal')
plt.show()
'''
