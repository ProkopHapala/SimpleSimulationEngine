#!/usr/bin/python
# use like:
# python replace.py answer.bas 1 14 2.0 0.75


import sys
import numpy as np
import atomicUtils as au

# ---- setup
fin       =     sys.argv[1]
typ       = int(sys.argv[2])
ofTyp     = int(sys.argv[3])
rcut      = float(sys.argv[4])
prob      = float(sys.argv[5])

reps = [
(9,1.60),
(17,2.02),
(35,2.15),
(53,2.43)
]

# bond lengths
#       E[kJ/mol]  Bond Length
# Si-F	565	       1.60
# Si-Cl	381   	   2.02
# Si-Br	310   	   2.15
# Si-I	234	       2.43
       
# ---- main

atoms   = np.genfromtxt( fin, skip_header=1 )

mask_Si = atoms[:,0]==ofTyp
mask_H  = atoms[:,0]==typ

for rep in reps:
    bond_counts = au/countTypeBonds( atoms, atoms[mask_H], rcut );
    mask_SiH1 = np.logical_and( ( bond_counts==1 ), mask_Si )
    found, foundDict = au.findBondsTo( atoms, 1, atoms[mask_SiH1], rcut=rcut )
    atoms2 = au.replace( atoms.copy(), found, to=rep[0], bond_length=rep[1], prob=prob  )
    au.saveAtoms( atoms2, "replaced_%i.xyz" %rep[0], xyz=True )
    au.saveAtoms( atoms2, "replaced_%i.bas" %rep[0], xyz=False )


