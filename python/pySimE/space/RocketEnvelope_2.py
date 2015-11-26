#!/usr/bin/python


def vCielkovsky( ve, mass_ratio ):
	return ve * log( 1/mass_ratio )

def sCielkovsky( t, dmdt, ve, velocity ):
	return ( t - 1/dmdt )*velocity + t*ve

def RocketEqTime( t, ve, specific_power, engine_ratio = 0.1, payload_ratio = 0.1, m0=1.0 ):
	'''
	t				[s]
	ve				[ m/s  ]
	payload_ratio	[ 1	  ]
	specific_power	[ W/kg ]
	'''
	construction_ratio =  ( payload_ratio + engine_ratio )
	specific_thrust = specific_power / ve                    # [ N / kg    = W / (m/s) / kg  ]
	thrust          = m0 * engine_ratio * specific_thrust    # N
	mass_flow       = thrust / ve                            # [ kg/s = N / (m/s)  ]
	dmdt = mass_flow/m0
	t_off = ( 1 - construction_ratio ) / dmdt


	mask = ( t > t_off )
	mass_ratio           = 1 - dmdt*t;

	acceleration         =  thrust / ( m0 * mass_ratio )
	acceleration[ mask ] = 0

	velocity_off         =  vCielkovsky( ve, construction_ratio )
	velocity             = 	vCielkovsky( ve, mass_ratio )
	velocity[ mask ]     = 	velocity_off 

	distance_off         =  sCielkovsky( t_off, dmdt, ve, velocity_off )
	distance             =  sCielkovsky( t    , dmdt, ve, velocity     )
	distance[ mask ]     =  distance_off + velocity_off*(t[ mask ] -t_off)

	print " thrust mass_flow dmdt t_off velocity_off distance_off : ", thrust, mass_flow, dmdt, t_off, velocity_off, distance_off
	return velocity, acceleration, distance

if __name__ == "__main__":
	from pylab import *
	figure(figsize=(15,5))
	t   = 10**( arange(-1,12,0.1) )
	v,a,s = RocketEqTime( t,    4462 , 9.81*136.7*4462 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);	# H2/LOX
	#v,a = RocketEqTime( t,    4462 , 9.81*136.7*4462, engine_ratio = 0.1, payload_ratio = 0.5 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);	# H2/LOX
	#v,a = RocketEqTime( t,    4462 , 9.81*136.7*4462, engine_ratio = 0.01, payload_ratio = 0.1 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s); 	# H2/LOX
	v,a,s = RocketEqTime( t,   8338.5,       10.0*8338.5 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);  # NERVA
	v,a,s = RocketEqTime( t, 189333.0  ,      250000/5.0 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);  # Ion_thruster DS4G ( Xe )
	v,a,s = RocketEqTime( t, 5170.0e+3  ,  1e+9/113400.0 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);  #  Werka FFRE
	show()
