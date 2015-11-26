# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 08:35:21 2013

@author: asiJa
"""
from pylab import *
from scipy.integrate import odeint



m0 = 1.0
#def func(m, x):
#	return m
	
#def func(m, x):
#	return m *( 10.0/x**2 + 1.0/x)
	
	
G      = 6.67384e-11  # N/kg^2/m^2
M      = 5.9736e+24   # kg
s      = 1988e+3      # N*m/kg  S-Glass
s      = 1790e+3      # N*m/kg  Basalt fiber
s      = 2071e+3      # N*m/kg  Vectran 
s      = 2071e+3      # N*m/kg  Carbon
s      = 2514e+3      # N*m/kg  Kevlar
s      = 3711e+3      # N*m/kg  Dyneema
#s      = 3766e+3 *2   # N*m/kg  Zylon PBO
s     = 5843e+3     # N*m/kg  PBO theoretical
#s     = 46268e+3    # N*m/kg  carbon nanotube
r0     = 6378.1e+3+25 # m
h0     = 1000e+3      # m
v0     = sqrt( G*M / (r0+h0)  )     # m/s
vair   = 2.5e+3      # m/s  
vc     = v0 - vair 

xs = arange(0.0,h0+0.0001, h0/100.0)

print v0

def vce(x):
	return v0 * (x+r0)/(r0+h0)

def vct(x):
	return vc * (h0-x) / h0

def a_cent(v, r):
	return v**2/r
	
def a_grav(x):
	return M*G/x**2
	
def func(m, x, s=3766e+3):
	ve = vce(x)
	vt = vct(x)
	return (m/s) *( a_grav(x+r0) + a_cent( vt,(h0+x)) - a_cent( ve-vt ,x+r0) )

def funcSimple(m, x):
	vt = vct(x)
	return (m/s) *( a_cent( vt,(h0+x)) )

def funcSimpleG(m, x):
	vt = vct(x)
	return (m/s) *( a_grav(x+r0) )

'''
ves = vce(xs)
vts = vct(xs)

#plot(xs,ves, 'x-r')
#plot(xs,vts, 'x-g')
#plot(xs,ves-vts, 'x-b')

aG  = a_grav(xs+r0)
aCt = a_cent( vts,(h0+xs))
aCe = a_cent( ves-vts ,xs+r0)

plot(xs,aG, 'x-r')
plot(xs,aCt,'x-g')
plot(xs,aCe,'x-b')
plot(xs,aG+aCt-aCe,'x-m')
axhline(0)
'''

ms  = odeint( func, m0, xs)                      #  zylon
ms2  = odeint( func, m0, xs, args=(5843e+3,))    #  PBO theoretical
ms3  = odeint( func, m0, xs, args=(46268e+3/2,)) #  carbon nanotube
#ms2 = odeint( funcSimple, m0, xs)
#ms3 = odeint( funcSimpleG, m0, xs)
plot(xs/1000,ms,'.-')
plot(xs/1000,ms2,'.-')
plot(xs/1000,ms3,'.-')
plot(xs/1000,zeros(len(xs)),'-'); xlabel('Height [km]'); ylabel('rope width')

#plot(xs,exp(xs))

show()

