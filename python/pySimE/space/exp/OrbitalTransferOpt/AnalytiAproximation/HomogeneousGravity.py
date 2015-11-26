#!/usr/bin/python



from __future__ import division
from sympy import *
#from mpmath import *

'''
This is approximate solution of movement of body with constant thrust in homogeneous gravity field

It is usefull aproximation to real orbital problem in case, when the curvature of gravitational potential is negligieble at scale of traveled distance

'''

r = Function('r')   # radial direction
s = Function('s')   # tangential direction
t, Ts, Tr, vr0, vt0 = symbols('t Ts Tr vr0 vt0')
GM, Gr0, Gr1, Gr2, Gs0, Gs1, Gs2 = symbols('GM, Gr0, Gr1, Gr2, Gs0, Gs1, Gs2')

R0, R, Fric = symbols('R0 R Fric')
#global_assumptions.add( Assume(Gr0, Q.positive) )
#global_assumptions.add( Assume(Gr1, Q.positive) )
#global_assumptions.add( Assume(Tr , Q.positive) )

EqS = s(t).diff(t,t) - Ts    ;                    print dsolve( EqS , s(t))

print diff( 1/(R0+R)**2, R )

r0 = pi
A = 1.0/(r0**2)
B = -2.0/(r0**3)

EqR = r(t).diff(t,t) - Tr  - A - B*r(t);      print dsolve( EqR , r(t))


EqS = s(t).diff(t,t) - Ts - pi*s(t).diff(t);                    print dsolve( EqS , s(t))


'''
solution

r(t) == C1*cos(sqrt(2)*t/R0**(3/2)) + C2*sin(sqrt(2)*t/R0**(3/2)) + 0.5*R0**3*Tr + 0.5*pi

'''

