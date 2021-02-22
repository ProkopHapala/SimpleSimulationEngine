#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

# ==== Setup

alpha = -1.6

atypes={
'H' : (1.487, 0.0006808),
'C' : (1.908, 0.0037292),
'O' : (1.661, 0.0091063),
'N' : (1.780, 0.0073719),
'Cl': (1.948, 0.0200000),
'Xe': (2.181, 0.0243442)
}

# ==== Functions


def getLJ( r, R, eps ):
    return eps*( (R/r)**12 - 2*(R/r)**6 )

def getMorse( r, R, eps, alpha=alpha ):
    return eps*( np.exp(2*alpha*(r-R)) - 2*np.exp(alpha*(r-R)) )

def getMorseParts(r, R, eps, alpha=alpha):
    EP = eps*np.exp(  2*alpha*(r-R) )
    EL = eps*np.exp(    alpha*(r-R) )
    #print EP, EL
    return EP, EL

def combMorseParts( EP, EL, R, eps, alpha=alpha ):
    CP =     eps*np.exp(-2*alpha*R)
    CL =  -2*eps*np.exp(-  alpha*R)
    print R, eps, alpha, CP, CL
    return ( EP*CP + EL*CL )

# ==== Main

xs = np.linspace(2.0,8.0,300)

elem1 = 'Cl'
atyp1 = atypes[elem1]

EP,EL = getMorseParts(xs, atyp1[0], np.sqrt(atyp1[1]) )

for elem2,clr in [ ('H','r'), ('Xe','b')]:
    atyp2 = atypes[elem2]
    Rij = atyp1[0]+atyp2[0]
    Eij = np.sqrt(atyp1[1]*atyp2[1])
    plt.plot( xs, getLJ(xs, Rij, Eij ),    ls='--', c=clr, label=('%s-%s LJ' %(elem1,elem2)) )
    plt.plot( xs, getMorse(xs, Rij, Eij ), ls=':' , c=clr, label=('%s-%s' %(elem1,elem2)) )
    plt.plot( xs, combMorseParts( EP, EL, atyp2[0], np.sqrt(atyp2[1]) ), ls='-', c=clr, label=('%s-%s factorized' %(elem1,elem2)) )
    plt.axvline( Rij, ls='--', c='k')
    

vmax=0.025
plt.ylim(-vmax,vmax)
plt.grid()
plt.legend()
plt.show()
