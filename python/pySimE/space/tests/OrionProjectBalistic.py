# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 13:48:15 2013

@author: asiJa
"""
from pylab import *

ktTNT = 4.184e+12 # Joule

EnergyDensity_U235 = 83.140e+12
EnergyDensity_DT   = 576e+12

const_eV 				= 1.60217656535e-19 # J
const_atomicMass			= 1.660538921e-27 # kg

EnergyDensity_U235_2  = 202.5e+6*const_eV /  ( const_atomicMass*235 )
EnergyDensity_DT_2    = 17.6e+6*const_eV  /  ( const_atomicMass*5   )
EnergyDensity_LiD_2   = 22.4e+6*const_eV  /  ( const_atomicMass*8   )
EnergyDensity_BeH_2   = 8.7e+6*const_eV   /  ( const_atomicMass*12  )

print " Energy denisty U235 : ", EnergyDensity_U235, EnergyDensity_U235_2
print " Energy denisty DT   : ", EnergyDensity_DT, EnergyDensity_DT_2
print  EnergyDensity_U235_2 / EnergyDensity_U235
print  EnergyDensity_DT_2 / EnergyDensity_DT

BombEnergy = 250*ktTNT;                                                       print " Energy ", BombEnergy, '[ J ]'
FissionMass = BombEnergy / EnergyDensity_U235_2;  BombMass=3*FissionMass;     print " Fission Mass ", FissionMass, '[kg] BombMass',BombMass,' [kg] '

def kineticVelocity(E,m):
	return sqrt( 2*E/m )

def E2impulse( E, m ):
	return 2*E/kineticVelocity(E,m)

v_exhaust = kineticVelocity(BombEnergy,BombMass);            print " v_exhaust ", v_exhaust /1000.0, '[ km/s ]'    
impulse   = E2impulse( BombEnergy*0.5, BombMass*0.5 );       print " Impulse   ", BombEnergy/1000.0, '[ kg*km/s ]'

m_pusher = 100e+3   # [kg]
m_ship = 1000e+3   # [kg]
v_pusher = (impulse/m_pusher)/2.0;           print 'v_pusher ='   ,v_pusher /1000.0,' [km/s]' 
v_ship   = (impulse/m_ship);                 print 'deltav_ship =',  v_ship /1000.0,' [km/s]'

specific_strength={ 'steel':254, 'Titanium':288,'Scifer':706, 'S-Glass':1988e+3, 'Basalt':1790e+3, 'Vectran':2071e+3, 'Carbon':2071e+3, 'Kevlar':2514e+3, 'Dyneema':3711e+3, 'Zylon':3766e+3, 'PBOmax':5843e+3 , 'CNT': 46268e+3 }

#   Slign_Shot Damper
L_rope = 5000
T_half = L_rope*pi/v_pusher;                         print 'T_half ='   ,T_half,' [ s ]' 
F_rope = m_pusher*(v_pusher**2)/L_rope;              print 'F_rope ='   ,F_rope /1000.0,' [ kN ]'   
M_rope = L_rope*F_rope/specific_strength['Carbon'];  print 'M_rope ='   ,M_rope ,' [ kg ]' 
ShipAcceleration  =  F_rope/m_ship;                  print 'ShipAcceleration ='   ,ShipAcceleration /9.81,' g ' 

def SlignOrion(L_rope=5000, m_ship= 1000e+3, m_pusher=100e+3, Bomb_Energy=250*ktTNT, Bomb_Mass=None, TransferEfficiency=0.5 ):
	if(Bomb_Mass==None):
		Bomb_Mass = 10 + 3*( BombEnergy / EnergyDensity_U235_2 )
	impulse  = 	sqrt(2*Bomb_Energy*Bomb_Mass*TransferEfficiency**2)
	v_pusher = (impulse/m_pusher)/2.0
	F_rope   = m_pusher*(v_pusher**2)/L_rope;
	ShipAcceleration  =  F_rope/m_ship
	return ShipAcceleration, T_half, v_pusher

print SlignOrion(Bomb_Mass=37.7428297645)

'''	
def SlignOrion( ShipAcceleration=9.81, L_rope=5000, m_ship= 1000e+3, m_pusher=100e+3, impulse=1.046e+12  ):
	if(L_rope=None): 
		v_pusher = (impulse/m_pusher)/2.0
	F_rope   = m_pusher*(v_pusher**2)/L_rope;
	ShipAcceleration  =  F_rope/m_ship
	return ShipAcceleration
'''



 


