#!/usr/bin/python

import constants as const
import numpy     as np

# ===========  Mechanics  ===================

def kineticVelocity( mass=const.protonMass, energy=1e+6*const.eV, relativistic=True  ):
	if relativistic:
		E0 = mass*const.lightSpeed**2
		return const.lightSpeed * pylab.sqrt( 1-(E0/(E0+energy))**2 ) # relativistic http://physics.stackexchange.com/questions/716/relativistic-speed-energy-relation-is-this-correct
	else:
		return pylab.sqrt( 2*energy/mass )
		
def kineticEnergy( mass=const.protonMass, velocity=0.1*const.lightSpeed, relativistic=True  ):
	if relativistic:
		E0 = mass*const.lightSpeed**2
		#return E0 / pylab.sqrt( 1.0 - (velocity/const.lightSpeed)**2 ) - E0
		return E0 * STRGamma(velocity) - E0
	else:
		return 0.5*mass*velocity**2


def classical_Energy_to_Velocity( Ek, m0 = const.protonMass ):
	return np.sqrt( 2*Ek/m0 )

def relativistic_Energy_to_Velocity( Ek, m0 = const.protonMass ):
	E0 = m0 * const.lightSpeed**2
	return const.lightSpeed * np.sqrt( Ek*(2*E0 + Ek )/(E0+Ek)**2 )

# ===========  Thermodynamics  ===================
	
def thermalRadiativePower( temperature = 3500, area = 1 ):
	return const.StefanBoltzmann*area*temperature**4
	
def equilibriumTemperature( power_density, absorptivity=1.0, emissivity=1.0 ):
	return ( power_density*(absorptivity/(const.StefanBoltzmann*emissivity) ) )**0.25

def heatCapacityRatio(degresOfFreedom=3):
	return 1 + (2.0/degresOfFreedom)

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

# ================ Electromagnetism =====================
	
def cyclotronRadius( velocity=0.05*const.lightSpeed, magneticB=1.0, mass = const.protonMass, charge=const.elementaryCharge, relativistic=True ):
	if relativistic:
		return STRGamma(velocity)*velocity*mass/(charge*magneticB)
	else:
		return velocity*mass/(charge*magneticB)

def cyclotronFrequency( magneticB=1.0, mass = const.protonMass, charge=const.elementaryCharge, velocity=None  ):
	if (velocity!=None):	# relativistic
		return charge*magneticB/( 2.0*pylab.pi*mass*STRGamma(velocity) )
	else:
		return charge*magneticB/( 2.0*pylab.pi*mass )

def HelmholzCoilField( current=1, radius=1 ):
	return 0.7155417528*cost_VacumPermeability*current/radius

# ============== Optics =====================

def difractionAperture( radius = 1.0, distance = 384.4e+6, waveLength = 1e-6 ):
	return 1.22*distance*waveLength/(2*radius)



