#!/usr/bin/python

import numpy as np

def mCielkovsky( ve, delta_v ):
	return np.exp( -delta_v/ve )

def vCielkovsky( ve, mass_ratio ):
	return ve * np.log( 1/mass_ratio )

def sCielkovsky( t, dmdt, ve, velocity ):
	return ( t - 1/dmdt )*velocity + t*ve

def payloadBoost( vexh, vend, vboost, mConstr=0.1 ):
    m_full         = mCielkovsky( vexh, vend )
    if (m_full < mConstr): 
        return ( vexh + vend + vboost ) * np.NAN
    m_boosted      = mCielkovsky( vexh, vend - vboost )
    return  ( m_boosted - mConstr ) / (m_full - mConstr)
    
def oberthHyperbolic( v1, vesc, delta_v ):
    e1    = (v1**2 + vesc**2)
    v1in  = np.sqrt(e1)
    v2in  = v1in + delta_v
    e2out = v2in**2 - vesc**2
    v2out = np.sqrt(e2out)
    return v2out
    
'''
def boostedMass( v0, vend, vexh ):
    #f_booster = mCielkovsky( vexh, v0       -vend*0 )
    f_full    = mCielkovsky( vexh, vend     -  v0*0 )
    f_boosted = mCielkovsky( vexh, vend-v0 )
    #print f_full, f_boosted, f_boosted * f_booster
    return f_full, f_boosted 
'''

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


# =================== Testing

def test_boostedMass():
    vs = np.linspace( 0.0,2000.0, 10)
    #mH2, mH2_   = boostedMass( vs, 7800.0, 4500.0 )
    #mCH4,mCH4_  = boostedMass( vs, 7800.0, 3615.0 )
    mH2LEO   = mCielkovsky( 4500.0, 7900.0  - vs )
    mH2escp  = mCielkovsky( 4500.0, 11186.0 - vs )
    mCH4LEO  = mCielkovsky( 3615.0, 7900.0  - vs )
    mCH4escp = mCielkovsky( 3615.0, 11186.0 - vs )
    
    mH2LEO   = payloadBoost( 4500.0, 7900.0,  vs, mConstr=0.1 )
    mH2escp  = payloadBoost( 4500.0, 11186.0, vs, mConstr=0.1 )
    mCH4LEO  = payloadBoost( 3615.0, 7900.0,  vs, mConstr=0.1 )
    mCH4escp = payloadBoost( 3615.0, 11186.0, vs, mConstr=0.1 )
    
    plt.plot(vs,mH2LEO   ,label='H2  LEO. ')
    plt.plot(vs,mH2escp  ,label='H2  escp.')
    plt.plot(vs,mCH4LEO  ,label='CH4 LEO. ')
    plt.plot(vs,mCH4escp ,label='CH4 escp.')
    mConstr=0.1
    #plt.plot(vs,mH2LEO-mConstr,label='H2  LEO. ')
    #plt.plot(vs,mH2escp  ,label='H2  escp.')
    #plt.plot(vs,mCH4LEO  ,label='CH4 LEO. ')
    #plt.plot(vs,mCH4escp ,label='CH4 escp.')
    plt.grid()
    plt.legend(loc=2)
    

def test_RocketEqTime():
    import matplotlib.pyplot as plt
    plt.figure(figsize=(15,5))
    t   = 10**( np.arange(-1,12,0.1) )
    v,a,s = RocketEqTime( t,    4462 , 9.81*136.7*4462 ); plt.subplot(1,3,1); plt.loglog(t,v); plt.subplot(1,3,2); plt.loglog(t,a); plt.subplot(1,3,3); plt.loglog(t,s);	# H2/LOX
    #v,a = RocketEqTime( t,    4462 , 9.81*136.7*4462, engine_ratio = 0.1, payload_ratio = 0.5 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s);	# H2/LOX
    #v,a = RocketEqTime( t,    4462 , 9.81*136.7*4462, engine_ratio = 0.01, payload_ratio = 0.1 ); subplot(1,3,1); loglog(t,v); subplot(1,3,2); loglog(t,a); subplot(1,3,3); loglog(t,s); 	# H2/LOX
    v,a,s = RocketEqTime( t,   8338.5,       10.0*8338.5 ); plt.subplot(1,3,1); plt.loglog(t,v); plt.subplot(1,3,2); plt.loglog(t,a); plt.subplot(1,3,3); plt.loglog(t,s);  # NERVA
    v,a,s = RocketEqTime( t, 189333.0  ,      250000/5.0 ); plt.subplot(1,3,1); plt.loglog(t,v); plt.subplot(1,3,2); plt.loglog(t,a); plt.subplot(1,3,3); plt.loglog(t,s);  # Ion_thruster DS4G ( Xe )
    v,a,s = RocketEqTime( t, 5170.0e+3  ,  1e+9/113400.0 ); plt.subplot(1,3,1); plt.loglog(t,v); plt.subplot(1,3,2); plt.loglog(t,a); plt.subplot(1,3,3); plt.loglog(t,s);  #  Werka FFRE
    plt.show()

def test_Oberth():
    import matplotlib.pyplot as plt
    dvs = [1.0e+3,2.0e+3,4.0e+3,8.0e+3]
    v1s = np.linspace( 0.0, 100.0e+3, 100 )
    plt.figure(figsize=(5.0,10.0))
    for dv in dvs:
        v2s = oberthHyperbolic( v1s, 60.0e+3, dv )
        plt.subplot(2,1,1); plt.plot( v1s*1e-3, v2s*1e-3,     label=("%3.1fkm/s"%(dv*1e-3))); 
        plt.subplot(2,1,2); plt.plot( v1s*1e-3, (v2s-v1s)/dv, label=("%3.1fkm/s"%(dv*1e-3)));
    plt.subplot(2,1,1); plt.grid(); plt.ylim(0,100.0); plt.axis('equal'); plt.xlabel("v1 [km/s]"); plt.ylabel("v2 [km/s]");
    plt.subplot(2,1,2); plt.grid(); plt.ylim(0,10.0);  plt.legend();      plt.xlabel("v1 [km/s]"); plt.ylabel("p=(v2-v1)/dv");
    plt.show()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #test_RocketEqTime()
    #test_boostedMass()
    test_Oberth()
    plt.show()
    

