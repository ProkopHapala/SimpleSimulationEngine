# -*- coding: utf-8 -*-
"""


How to debug:
	from debug import *
	import KosmoSuite as ks
	ks.equilibriumTemperature(1000)
	... modify something in equilibriumTemperature()
	reset_sys(); import KosmoSuite as ks
	ks.equilibriumTemperature(1000)


example of use:
    
>>> import KosmoSuit_main as ks
>>> ks.exhaustVelocity(molarMass=18, kappa=1.2222, temperature=3500)

"""
from __future__ import division

import pylab

const_StefanBoltzmann        = 5.67037321e-8     # W m^-2 K^−4
const_universalGas           = 8.3144621         # J /K /mol 
const_atmosphericPressure    = 101325.0          # Pa
const_gravitational          = 6.67384e-11       # m^3 kg / s^2
const_electronCharge         = 1.60217657e-19    # coulombs
const_electronMass           = 9.109382910e-31   # kg
const_atomicMass             = 1.660538921e-27   # kg
const_protonMass             = 1.672621780e-27   # kg
const_neutronMass            = 1.674927351e-27   # kg
const_lightSpeed             = 299792458.0	     # m/s
const_AU                     = 149597871.0e+3    # m
const_0celsius               = 273.15            # K   
cost_VacumPermeability       = 4e-7*pylab.pi 	 # Tesla·m/A
const_eV                     = 1.60217656535e-19 # J

#execfile( 'table_ChemicalFuels.py' )


# ========== Electrical ====================

def conductorMass( length, power=1, voltage=1, loss=0.5, resistivity=2.82e-8, density=2.70e+3 ):
	# P_tot = P_loss + P_used
	# U_tot = U_loss + U_used | . I
	# R_tot = R_loss + R_used | . I^2
	#I = power/voltage
	#R_loss = (loss*power)/I**2 
	#S = resistivity*length/R_loss
	S = resistivity*length*power/(loss*voltage**2)
	#print I,R_loss,S
	return S * length * density
	
# ========== Chemistry =====================

def oxygenShare(CO2=0,H2O=0,CO=0):
	fuel 		= 12*(CO2+CO) + 2 *H2O
	exhaust	= 44*(CO2+CO) + 18*H2O
	return fuel/exhaust

# ===========  Thermodynamics  ===================
	
def thermalRadiativePower( temperature = 3500, area = 1 ):
	return const_StefanBoltzmann*area*temperature**4
	
def equilibriumTemperature( power_density, absorptivity=1.0, emissivity=1.0 ):
	return ( power_density*(absorptivity/(const_StefanBoltzmann*emissivity) ) )**0.25

def heatCapacityRatio(degresOfFreedom=3):
	return 1 + (2.0/degresOfFreedom)

# ==========

# Example for H2/O2:        exhaustVelocity(molarMass=18, kappa=1.2222, temperature=3500)   = 4217.2462317444688
# Example for H2 nuclear	exhaustVelocity(molarMass=2, kappa=1.4, temperature=3500)   = 10092.183149596525
def exhaustVelocity( energyDensity=None, pressureOut=None, pressureIn=const_atmosphericPressure*25, temperature=3500, kappa=1.666666, molarMass=2 ):
	if (energyDensity!=None):	# from kinetic energy
		return pylab.sqrt( 2*energyDensity )
	else:	 # lavalNozzle adiabatic expansion http://en.wikipedia.org/wiki/De_Laval_nozzle
		gamma = (kappa)/(kappa-1)
		ExpansionTerm = 1.0
		if (( pressureOut != None ) ):
				ExpansionTerm = (1.0 - (pressureOut/pressureIn)**(1/gamma) )
		return pylab.sqrt(2000.0*const_universalGas*temperature*gamma*ExpansionTerm/molarMass)			

def tsilkovskyPayload( deltaV=9400, exhaustVelocity=4462 ):
	return pylab.exp(-float(deltaV)/exhaustVelocity)
	
def tsilkovskyVelocity( massRatio=0.1, exhaustVelocity=4462 ):
	return exhaustVelocity*-pylab.log(massRatio)

def deltaVOberth( deltaV = 10, escapeVeloctiy = 617.5 ):
	return pylab.sqrt(deltaV**2 + 2.0*escapeVeloctiy*deltaV)	

# ================== Orbital Dynamics ==============

def gravitationalAcceleration( mass=5.9736e+24 , radius=6378100 ):
	return const_gravitational*mass/float(radius)**2

def orbitalVelocity(  mass=5.9736e+24 , radius=6378100 ):
	return pylab.sqrt( const_gravitational*mass / float(radius) )	

# example:    escapeVelocit = sqrt(2)orbitalVelocity
# example:    escapeVelocity() = 11180.862334866215	for earth from earth surface	 	
def escapeVelocity( mass=5.9736e+24 , radius=6378100, denisty=None ):
	if (denisty!=None):
		return 2.364e-5 * float(radius)*pylab.sqrt(density)   # for spherical body
	else:	
		return pylab.sqrt( 2.0*const_gravitational*mass / radius  )	

# example:    GEO  orbitalPeriod( semimajorAxis=42e+6  ) / 24  = 	3568.9182312398239
def orbitalPeriod( semimajorAxis = 6378100, mass=5.9736e+24 ):
	return pylab.pi*2.0*pylab.sqrt( semimajorAxis**3.0 / ( const_gravitational*mass ) )
	
def centrifugalAcceleration( velocity=7906, radius=6378100, omega=None):
	if (omega!=None):
		return (omega**2)*radius
	else:	
		return velocity**2/radius

STRGamma = lambda v: 1.0/pylab.sqrt( 1.0-(v/const_lightSpeed)**2)  

def kineticVelocity( mass=const_protonMass, energy=1e+6*const_eV, relativistic=True  ):
	if relativistic:
		E0 = mass*const_lightSpeed**2
		return const_lightSpeed * pylab.sqrt( 1-(E0/(E0+energy))**2 ) # relativistic http://physics.stackexchange.com/questions/716/relativistic-speed-energy-relation-is-this-correct
	else:
		return pylab.sqrt( 2*energy/mass )
		
def kineticEnergy( mass=const_protonMass, velocity=0.1*const_lightSpeed, relativistic=True  ):
	if relativistic:
		E0 = mass*const_lightSpeed**2
		#return E0 / pylab.sqrt( 1.0 - (velocity/const_lightSpeed)**2 ) - E0
		return E0 * STRGamma(velocity) - E0
	else:
		return 0.5*mass*velocity**2

# ================ Electromagnetism =====================
	
def cyclotronRadius( velocity=0.05*const_lightSpeed, magneticB=1.0, mass = const_protonMass, charge=const_electronCharge, relativistic=True ):
	if relativistic:
		return STRGamma(velocity)*velocity*mass/(charge*magneticB)
	else:
		return velocity*mass/(charge*magneticB)

def cyclotronFrequency( magneticB=1.0, mass = const_protonMass, charge=const_electronCharge, velocity=None  ):
	if (velocity!=None):	# relativistic
		return charge*magneticB/( 2.0*pylab.pi*mass*STRGamma(velocity) )
	else:
		return charge*magneticB/( 2.0*pylab.pi*mass )

def HelmholzCoilField( current=1, radius=1 ):
	return 0.7155417528*cost_VacumPermeability*current/radius

# ============== plazma ======================

def lossCone( fieldRatio = 2.0 ):
	return 1.0/pylab.arcsin( pylab.sqrt( fieldRatio ) )


# ============== Optics =====================

def difractionAperture( radius = 1.0, distance = 384.4e+6, waveLength = 1e-6 ):
	return 1.22*distance*waveLength/(2*radius)

# ============== Space Elevator ============

def spaceElevatorS( sigma=64e+9, rho=1340.0, omega = 7.2921150e-5, g0=9.780, r0=6371000.0  ):
	x = ( omega**2 )*r0*g0
	return pylab.exp( (rho/sigma) *g0*r0*(1+0.5*x-1.5*(x**0.3333333)))	
	
# ================ RailGun =================

