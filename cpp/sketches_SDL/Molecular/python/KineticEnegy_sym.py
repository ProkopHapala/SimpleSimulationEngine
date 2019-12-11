#!/usr/bin/python

import sympy as sy

'''

Increase of Kinetic energy by orthogonalization of two everlaping orbitas f1,f2

d1,d2 are laplacians of the orbitals

we have two linear combinations
    fsyn  = f1 + f2
    fanti = f1 - f2
and their laplacians
    dsyn  = d1 + d2
    danti = d1 - d2

This leads to matrix elements of kinteic energy re-normalized to make sure that combined wavefunction has 1 electron each
    Tsyn  = fsyn*dsyn/(fsyn*fsyn)
    Tanti = fanti*danti/(fanti*fanti)

This can be efficiently expressed using bare matrix elements of 
kinetic energy 
    t11 = f1*d1
    t12 = f1*d2
    t21 = f2*d1
    t22 = f2*d2
and overlap:
    s11 = 1
    s22 = 1
    s = s12 = s21 = f2*f1 = f1*f2

Tsyn  = (t11+t22+t12+t21)/(s11+s22+s12+s21) = (t11+t22+t12+t21)/(2(1+s))
Tanti = (t11+t22-t12-t21)/(s11+s22-s12-s21) = (t11+t22-t12-t21)/(2(1-s))

Total kinetic energy difference due to orthogonalization is:

dT = Tsyn + Tanti - t11 - t22
dT = s*( s*(t11 + t22) - t12 - t21)/( 1 - s**2 )


dT  = (s^2)*(E1 + E2)/( 1 - s^2 )    -    s*(t12 + t21)/( 1 - s^2 )

Assuming s^2 is very small

dT  ~= s*(t12 + t21)/( 1 - s^2 )


'''



'''
f1, f2, d1, d2 = sy.symbols('f1 f2 d1 d2')

fsyn  = f1 + f2
fanti = f1 - f2
dsyn  = d1 + d2
danti = d1 - d2

t11 = f1*d1
t12 = f1*d2
t21 = f2*d1
t22 = f2*d2

#s11 = f1*f1
#s22 = f2*f2
s11 = 1
s22 = 1
s12 = f1*f2
s21 = f2*f1
'''

t11, t22, t12, t21, s = sy.symbols('t11 t22 t12 t21 s')

s11 = 1
s22 = 1
s12 = s
s21 = s

Tsyn  = (t11+t22+t12+t21)/(s11+s22+s12+s21)
Tanti = (t11+t22-t12-t21)/(s11+s22-s12-s21)

# ========= Check  Tsyn  Tanti
#Tsyn_  = fsyn*dsyn/(fsyn*fsyn)
#Tanti_ = fanti*danti/(fanti*fanti)
#print "check Tsyn:  ", sy.simplify( Tsyn_  - Tsyn   )
#print "check Tanti: ", sy.simplify( Tanti_ - Tanti  )

dT = Tsyn + Tanti - t11 - t22

dT_ = s*(s*(t11 + t22) - t12 - t21)/( 1 - s**2 )

print "change of kinetic energy: "
print "DeltaT ", sy.simplify( dT )
print "DeltaT ", sy.simplify( dT_ )
print "check : ", sy.simplify( dT_ - dT )