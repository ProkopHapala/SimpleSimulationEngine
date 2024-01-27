import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from   matplotlib import collections  as mc

#os.system('mode con: cols=100 lines=50')
np.set_printoptions(linewidth=np.inf)

# =========== Functions

def drawStricks( ps, sticks, ax=None, color='k', lw=1, ls='-' ):
    if ax is None: ax = plt.gca()
    lines = []
    for i1,i2 in sticks:
        p1 = ps[i1][:2]
        p2 = ps[i2][:2]
        lines.append( [p1,p2] )
    lc = mc.LineCollection(lines, colors=color, linewidths=lw, linestyle=ls  )
    ax.add_collection(lc)
    #ax.axis('equal')

def build_Jacobian( sticks, ps ):
    nst = len(sticks)
    nps = len(ps)
    l0s = np.zeros( len(sticks) )
    hs  = np.zeros( (len(sticks),3) )  
    for ib,(i,j) in enumerate(sticks):
        d       = ps[i] - ps[j]
        l       = np.linalg.norm( d )
        l0s[ib] = l
        hs [ib] = d/l
    J = np.zeros( (nst, 3*nps) )
    for ib,(i,j) in enumerate(sticks):
        J[ib,3*i:3*i+3] = -hs[ib]
        J[ib,3*j:3*j+3] = +hs[ib]
    return J, l0s, hs

def build_MassMatrix( ms ):
    nps = len(ms)
    M_diag = np.zeros(3*nps)
    for i, m in enumerate(ms):
        M_diag[i*3:i*3+3] = m
    M = np.diag(M_diag)
    return M

# =========== Setup the system
    
# Masses of the points
ms = np.array([
10.0,  # Mass of point 1 (heavier)
1.0 ,  # Mass of point 2
1.5 ,  # Mass of point 3
])

# Positions of the points
ps = np.array([ 
    [ 0.0,  0.0, 0.0],
    [+1.0, -1.0, 0.0],
    [-1.0, -1.0, 0.0],
])

# Sticks connecting the points
sticks = [
    (0,1),
    (1,2),
    (2,0),
]

# =========== Main Body

J, l0s, hs = build_Jacobian  ( sticks, ps )
M          = build_MassMatrix( ms )

# Compute the matrix A  system
A = J @ M @ J.T

# Output the final matrix J and A to the console
print("Matrix J:\n", J)
print("Matrix A:\n", A)

plt.figure(figsize=(10,10))
drawStricks( ps, sticks, ax=None, color='k', lw=1, ls='-' )
plt.plot( ps[:,0], ps[:,1], 'o')
plt.grid()
plt.axis('equal')
plt.show()