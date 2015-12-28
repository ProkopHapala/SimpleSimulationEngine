#!/usr/bin/env python

from sympy import *
from sympy.matrices import Matrix

print 
print " ==== Solution of 6th order polygon coeficients ==== "

p0, p1, v0, v1, a0,a1,t = symbols('p0 p1 v0 v1 a0 a1 t')
"      1   t  t2  t3  t4  t5   "
A = Matrix([
    [  1,  0,  0,  0,  0,  0  ],  # p0
    [  1,  1,  1,  1,  1,  1  ],  # p1
    [  0,  1,  0,  0,  0,  0  ],  # v0
    [  0,  1,  2,  3,  4,  5  ],  # v1
    [  0,  0,  2,  0,  0,  0  ],  # a0
    [  0,  0,  2,  6,  12, 20 ]]) # a1
b = Matrix(6,1,[p0,p1,v0,v1,a0,a1])
CP = A.LUsolve(b)
print "Coefficients "
print CP

print 
print " ==== Make A4 A5 zero ==== "
A = Matrix([
    [  1.5,-1.0   ],
    [ -0.5, 0.5 ]])

b = Matrix(2,1,[
    -(15*p0 - 15*p1 + 8*v0 + 7*v1),
    -(-6*p0 +  6*p1 - 3*v0 - 3*v1)])

c = A.LUsolve(b)
print " Coefficients: "
print c

print 
print " ==== Make Derivative of accelaration zero === "
A3 = -1.5*a0 + 0.5*a1 - 10*p0 + 10*p1 - 6*v0 - 4*v1
A4 =  1.5*a0 -     a1 + 15*p0 - 15*p1 + 8*v0 + 7*v1
A5 = -0.5*a0 + 0.5*a1 -  6*p0  + 6*p1 - 3*v0 - 3*v1

#Eq1 = 6*A3
#Eq2 = 6*A3 + 24*A4 + 60*A5

Eq1 = 6*CP[3]
Eq2 = 6*CP[3] + 24*CP[4] + 60*CP[5]

aa = solve( [Eq1,Eq2],[a0,a1]  )
print aa


print 
print " ==== Polynom itselfs === "
P   = CP[0] + CP[1]*t +   CP[2]*t**2 +   CP[3]*t**3 +    CP[4]*t**4 +    CP[5]*t**5 
dP  =         CP[1]*t + 2*CP[2]*t    + 3*CP[3]*t**2 +  4*CP[4]*t**3 +  5*CP[5]*t**4 
ddP =                   2*CP[2]      + 6*CP[3]*t    + 12*CP[4]*t**2 + 20*CP[5]*t**3 
print  "P   =  ",P
print
print  "dP  =  ",dP
print 
print  "ddP =  ",ddP

print 
print " ==== Variational derivatives of gravitational acclearation === "
print " ::    F       =    1/r^2   ;   r = P(t)  "
print " :: dF(t)/dCi  =   -2/(P(t)^3) * dP/dCi   "

dP_da0 = diff( P, a0, 1)
print  "dP/da0 = ", dP_da0 
dP_da1 = diff( P, a1, 1)
print  "dP/da1 = ", dP_da1

print 
print " ==== Variational derivatives of kinematic acclearation in radial coordinates === "
# http://en.wikipedia.org/wiki/Polar_coordinate_system  Vector calculus
#   X =   
#  dX =                   dR*uR   +              R*dPhi*uPhi
# ddX = ( ddR + R*dPhi**2  )*uR   + R*ddPhi   + 2* dR * dPhi

print "  F_Phi = R * ddPhi     + 2 * dR * dPhi   "
print "  F_R   = R *  dPhi**2  +    ddR          "

ap0, ap1, ar0, ar1 = symbols('ap0, ap1, ar0, ar1')

Phi   = Function('Phi'  )(ap0,ap1)
dPhi  = Function('dPhi' )(ap0,ap1)
ddPhi = Function('ddPhi')(ap0,ap1)

R   = Function('R'  )(ar0,ar1)
dR  = Function('dR' )(ar0,ar1)
ddR = Function('ddR')(ar0,ar1)

F_Phi = R * ddPhi     + 2 * dR * dPhi
F_R   = R *  dPhi**2  +    ddR 

print
print " F_Phi = ", F_Phi

print 
print " d F_Phi / ap0 =  "
#pprint ( diff(  F_Phi, ap0) )
print  diff(  F_Phi, ap0) 
print " d F_Phi / ar0 =  "
#pprint ( diff(  F_Phi, ar0))
print  diff(  F_Phi, ar0)

print
print " F_Phi = ", F_Phi

print
print " d F_R / ap0   =  "
#pprint ( diff(  F_R, ap0) )
print  diff(  F_R, ap0) 
print " d F_R / ar0   =  "
#pprint ( diff(  F_R, ar0) )
print  diff(  F_R, ar0) 

diff_P_a0   = diff(   P, a0)
diff_P_a1   = diff(   P, a1)

diff_dP_a0  = diff(  dP, a0)
diff_dP_a1  = diff(  dP, a1)

diff_ddP_a0 = diff( ddP, a0)
diff_ddP_a1 = diff( ddP, a1)

print
print " diff_P_a0     =  ",diff_P_a0
print
print " diff_P_a1     =  ",diff_P_a1

print
print " diff_dP_a0    =  ",diff_dP_a0
print
print " diff_dP_a1    =  ",diff_dP_a1

print
print " diff_ddP_a0   =  ",diff_ddP_a0
print
print " diff_ddP_a1   =  ",diff_ddP_a1


'''
print 
print " ==== FULL EXAPNSION of kinematic variational derivatives === "

h0,h1,dh0,dh1,ddh0,ddh1 = symbols('h0 h1 dh0 dh1 ddh0 ddh1')
r0,r1,dr0,dr1,ddr0,ddr1 = symbols('r0 r1 dr0 dr1 ddr0 ddr1')

listP   = [ p0 ,p1 ,v0  ,v1  ,a0   ,a1   ]
listPhi = [ h0 ,h1 ,dh0 ,dh1 ,ddh0 ,ddh1 ] 
listR   = [ r0 ,r1 ,dr0 ,dr1 ,ddr0 ,ddr1 ]

PhiP = P.subs( listP[0], listPhi[0]  )
for i in range(1,len(listPhi)):
    PhiP = PhiP.subs( listP[i], listPhi[i] )

RP = P.subs( listP[0], listR[0]  )
for i in range(1,len(listR)):
    RP = RP.subs( listP[i], listR[i] )

dPhiP  = diff( PhiP, t, 1)
ddPhiP = diff( dPhiP, t, 1)
dRP    = diff( RP, t, 1)
ddRP   = diff( dRP, t, 1)
F_Phi  = RP * ddPhiP    +   2*dRP * dPhiP
F_R    = RP *  dPhiP**2 +    ddRP


print
print " ===== Derivatives of F_Phi  "

dFPhi_ddh0 = diff( F_Phi, ddh0, 1)
dFPhi_ddh1 = diff( F_Phi, ddh1, 1)
dFPhi_ddr0 = diff( F_Phi, ddr0, 1)
dFPhi_ddr1 = diff( F_Phi, ddr1, 1)

print
print " dFPhi_ddh0 =  ",dFPhi_ddh0
print
print " dFPhi_ddh1 =  ",dFPhi_ddh1
print
print " dFPhi_ddr0 =  ",dFPhi_ddr0
print
print " dFPhi_ddr1 =  ",dFPhi_ddr1

print
print " ===== Derivatives of F_R  "

dFR_ddh0   = diff( F_R  , ddh0, 1)
dFR_ddh1   = diff( F_R  , ddh1, 1)
dFR_ddr0   = diff( F_R  , ddr0, 1)
dFR_ddr1   = diff( F_R  , ddr1, 1)


print
print " dFR_ddh0 =  ",dFR_ddh0
print
print " dFR_ddh1 =  ",dFR_ddh1
print
print " dFR_ddr0 =  ",dFR_ddr0
print
print " dFR_ddr1 =  ",dFR_ddr1
'''
