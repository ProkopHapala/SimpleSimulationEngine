# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 19:47:14 2013

@author: asiJa
"""


# ==============================================

print 'space elevator'

from sympy import *
import pylab

GM,v,h0,r0,x, s, omega  = symbols('GM v h0 r0 x s omega')

m = Function('m')(x)


#print dsolve(Derivative(m, x) - m(x)*((GM/((r0+x)**2))+(v**2/(h0-x)))/s , m(x) )
#print dsolve(Derivative(m, x) - m(x)*(   (GM/(x**2)) + (v**2/(r0-x))   )/s , m(x) )
print dsolve(Derivative(m, x) - m(x)*(   (GM/(x**2)) - ((omega**2)*x)   )/s , m(x) )


G     = 6.67384e-11 # N/kg^2/m^2
M     = 5.9736e+24  # kg
s     = 3766e+3     # N*m/kg
r0    = 6378100+1000e+3     # m
v     = 6.89e+3     # m/s
omega = 7.2921150e-5 
x0 = 6378100


def constG(x):
	return pylab.exp( + 9.81*x/s)

def justGravity(x):
	return pylab.exp( - ( G*M/x ) / s)
	
def Full(x):
	print type(x)
	if (type(x) == pylab.ndarray):
		for xx in x:
			grav = G*M/xx
			cent = ((xx*omega)**2)/2
			print grav, cent, grav/cent
	return pylab.exp( - ( G*M/x + ((x*omega)**2)/2 ) / s)


xs = pylab.arange( 1,36000e+3,1000e+3 )





# Zylon
s     = 3766e+3     # N*m/kg
msConst = constG     (xs   )
msGrav  = justGravity(xs+x0)/justGravity(x0)
msFull  = Full       (xs+x0)/Full       (x0)
#pylab.plot(xs,msConst,'r-x')
pylab.plot(xs/1000,msGrav,'g-o')  #
pylab.plot(xs/1000,msFull,'g-x')  #

# Carbon nanotube
s     = 46268e+3     # N*m/kg
msConst = constG     (xs   )
msGrav  = justGravity(xs+x0)/justGravity(x0)
msFull  = Full       (xs+x0)/Full       (x0)
pylab.plot(xs/1000,msGrav,'r-o')  #
pylab.plot(xs/1000,msFull,'r-x')  #

pylab.yscale('log')
pylab.show()