#!/usr/bin/python

# compute distance by integration of tsielkovsky equation
# http://physics.stackexchange.com/questions/161358/how-can-i-caculate-the-time-to-traverse-a-distance-with-the-tsiolkovsky-equation
# v(t) =    v0        +  vexh log( mtot /( mtot - f*t ) ) 
# s(t) =  ( t - mtot/dm ) * vexh log( mtot /( mtot - dm*t ) )    +    t*(v0 + ve ) 
# s(t) =  tburn * ( t/tburn - mtot/mprop ) * vexh log( mtot /( mtot - dm*t ) )    +    t*(v0 + ve ) 


import numpy as np
import matplotlib.pylab as plt

# ========== Functions

def timeToDistace( times, specPower, vexh, mPayload, mEngine, mPropelant, nodes = np.linspace( 0.0,1.0,10 ) ):
	mTot     = mPayload + mEngine + mPropelant
	mDry     = mPayload + mEngine
	power    = specPower * mEngine 
	thrust   = power / vexh
	massFlow = thrust / vexh
	burnTime = mPropelant / massFlow
	ts       = nodes * burnTime
	vt       = vexh * np.log( mTot /( mTot - massFlow*ts ) )
	st       = vexh * ( ts - mTot/massFlow ) * np.log( mTot /( mTot - massFlow*ts ) )   +    ts*vexh
	if times[-1] > ts[-1]:
		st = np.append( st,    st[-1] + vt[-1]*(times[-1] - ts[-1])  )
		vt = np.append( vt,    vt[-1]                                )
		ts = np.append( ts,    times[-1]                             )
	st_ = np.interp( times,  ts, st  )
	vt_ = np.interp( times,  ts, vt  )
	mask     = times < ts[1]
	st_[mask] = 0.5 * (thrust/mTot) * times[mask]**2
	return st_, vt_

def combineSpecPower( thruster, powerPlant ):
	combined=  1.0/( 1.0/thruster  +   1.0/powerPlant )
	#print( "thruster, powerPlant, combined",   thruster, powerPlant, combined )
	return combined

# ========== Evaluation of ships

'''
=============== Power plants =====================


Jet engine			           10.0e+3  W/kg 			https://en.wikipedia.org/wiki/Power-to-weight_ratio        
Baryton turbine max            150.0e+3 W/kg			https://en.wikipedia.org/wiki/Power-to-weight_ratio
Brushless generator            10.0e+3  W/kg			https://en.wikipedia.org/wiki/Power-to-weight_ratio
High power density generator   10.0e+3  W/kg			DOI: 10.1109/IECEC.1989.74606


				Power[W]    Weight       W/kg 
NERVA
MITEE           0.3400E+6
TOPAZ			5.0e+3		320.0	    15.625					https://en.wikipedia.org/wiki/TOPAZ_nuclear_reactor
Bimodal 		25.0e+3		2224.0      11.24
phovolto        40.0		8.0e-3    5000.0

=> MITEE + Baryton turbine + High power density generator 
conservative     3e+3  W/kg
optimistic      30e+3  W/kg   

=============== Thrusters plants =====================

				Power[W]    Mass[W]   Exhaust
VASIMR_high     	6.0e+6	10000kg		294e+3           http://www.projectrho.com/public_html/rocket/enginelist.php#id--Electrothermal--Resistojet
VASIMR_med 			6.0e+6 	10000kg		147e+3
VASIMR_low 			6.0e+6 	10000kg		 29e+3
Mass driver         0.3e+12  150e+3 	 30e+3  
Coloide Truster     172e+6    20e+3		 43e+3
Ion thruster Cs    1050e+6   400e+3		210e+3
'''

nucConserv           = combineSpecPower(  10e+3  , combineSpecPower(  330e+3/5  ,  150e+3  ) )
nucOptimist          = combineSpecPower(  150e+3 , combineSpecPower(  330e+3/5  ,  150e+3  ) )

BimodalVASMIR        = combineSpecPower(   6.0e+6/10.0e+3 ,  25.0e+3/2224.0 )
photovoltVASIMR      = combineSpecPower(   6.0e+6/10.0e+3 ,     40.0/8.0e-3 )
photoMassDriver      = combineSpecPower(  0.3e+12/150e+3  ,     40.0/8.0e-3 )
photoColoid          = combineSpecPower(   172e+6/20e+3   ,     40.0/8.0e-3 )
photoCesium          = combineSpecPower(  1050e+6/400e+3  ,     40.0/8.0e-3 )

nucConservMassDriver = combineSpecPower(  0.3e+12/150e+3  ,     nucConserv )
nucConservColoid     = combineSpecPower(   172e+6/20e+3   ,     nucConserv )
nucConservCesium     = combineSpecPower(  1050e+6/400e+3  ,     nucConserv )

'''
ships=[
#   spec.Power [W/kg]    vexh [m/s]  payload   engine   propelant  description                            color      line 
[   1.2711e+6,         4400.0,     0.25,     0.25,    0.50,      'SSME  H2+LOX           4.4km/s 1:1: 2',  '#FF0000', '-' ],
[   0.3400e+6,         9800.0,     0.25,     0.25,    0.50,      'MITEE H2               9.8km/s 1:1: 2',  '#FF8000', '-' ],
[   photoMassDriver,      30e+3,     0.25,     0.25,    0.50,      'MassDriverPhovoltaic  30  km/s 1:1: 2',  '#8080FF', '-' ],
[   photoColoid,      	  43e+3,     0.25,     0.25,    0.50,      'ColoidElecPhovoltaic  43  km/s 1:1: 2',  '#008000', '-' ],
[   photoCesium,         210e+3,     0.25,     0.25,    0.50,      'Cs+IonElecPhovoltaic 200  km/s 1:1: 2',  '#0080FF', '-' ],
[   1.2711E+006,         4400.0,     0.05,     0.05,    0.90,      'SSME  H2+LOX           4.4km/s 1:1:18', '#FF0000', '--'  ],
[   0.3400E+006,         9800.0,     0.05,     0.05,    0.90,      'MITEE H2               9.8km/s 1:1:18', '#FF8000', '--'  ],
[   photoMassDriver,      30e+3,     0.05,     0.05,    0.90,      'MassDriverPhovoltaic  30  km/s 1:1:18', '#8080FF', '--'  ],
[   photoColoid,      	  43e+3,     0.05,     0.05,    0.90,      'ColoidElecPhovoltaic  43  km/s 1:1:18', '#008000', '--'  ],
[   photoCesium,         210e+3,     0.05,     0.05,    0.90,      'Cs+IonElecPhovoltaic 200  km/s 1:1:18', '#0080FF', '--'  ],
]
'''


'''
ships=[
#   spec.Power [W/kg]    vexh [m/s]  payload   engine   propelant  description                            color      line 
[   1.2711e+6,         		4400.0,     0.25,     0.25,    0.50,      'SSME  H2+LOX           4.4km/s 1:1: 2',  '#FF0000', '-' ],
[   0.3400e+6,         		9800.0,     0.25,     0.25,    0.50,      'MITEE H2               9.8km/s 1:1: 2',  '#FF8000', '-' ],
[   nucConservMassDriver,    30e+3,     0.25,     0.25,    0.50,      'MassDriverNucConserv  30  km/s 1:1: 2', '#8080FF', '-' ],
[   nucConservColoid,      	 43e+3,     0.25,     0.25,    0.50,      'ColoidElecNucConserv  43  km/s 1:1: 2', '#008000', '-' ],
[   nucConservCesium,       210e+3,     0.25,     0.25,    0.50,      'Cs+IonElecNucConserv 200  km/s 1:1: 2', '#0080FF', '-' ],
[   1.2711E+006,            4400.0,     0.05,     0.05,    0.90,      'SSME  H2+LOX           4.4km/s 1:1:18', '#FF0000', '--'  ],
[   0.3400E+006,          	9800.0,     0.05,     0.05,    0.90,      'MITEE H2               9.8km/s 1:1:18', '#FF8000', '--'  ],
[   nucConservMassDriver,    30e+3,     0.05,     0.05,    0.90,      'MassDriverNucConserv  30  km/s 1:1:18', '#8080FF', '--'  ],
[   nucConservColoid,      	 43e+3,     0.05,     0.05,    0.90,      'ColoidElecNucConserv  43  km/s 1:1:18', '#008000', '--'  ],
[   nucConservCesium,       210e+3,     0.05,     0.05,    0.90,      'Cs+IonElecNucConserv 200  km/s 1:1:18', '#0080FF', '--'  ],
]
'''

ships=[
#   spec.Power [W/kg]    vexh [m/s]  payload   engine   propelant  description                            color      line 
[   1.2711e+6,         		4400.0,     0.1,     0.1,    0.80,      'SSME  H2+LOX           4.4km/s 1:1:8',  '#FF0000', '-' ],
[   0.3400e+6,         		9800.0,     0.1,     0.1,    0.80,      'MITEE H2               9.8km/s 1:1:8',  '#FF8000', '-' ],
[   nucConservMassDriver,    30e+3,     0.1,     0.1,    0.80,      'MassDriverNucConserv  30  km/s 1:1:8', '#8080FF', '-' ],
[   nucConservColoid,      	 43e+3,     0.1,     0.1,    0.80,      'ColoidElecNucConserv  43  km/s 1:1:8', '#008000', '-' ],
[   nucConservCesium,       210e+3,     0.1,     0.1,    0.80,      'Cs+IonElecNucConserv 200  km/s 1:1:8', '#0080FF', '-' ],
]

#  =========== Plot formating

timeLines = np.array([1.0,  60.0, 3600.0, 86400.0, 604800,  2592000, 31556926, 315569260,  3155692600, 31556926000 ])
timeTexts          = ['sec','min','hour', 'day',   'week',  'month', 'year',   '10years',  '100years', '1000years']

#distLines = np.array([1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 6371e+3, 42164e+3, 384400e+3, 1e+9, 1e+10, 5.790918E+010, 1.082089E+011, 1.495979E+011, 2.279366E+011, 7.784120E+011, 1.426725E+012, 2.870972E+012,  4.498253E+012, 1.40621998e+13, 2.99195741e+14, 7.47989354e+15, 4.13425091e+16])
#distTexts =          ['10m','100m', '1km','10km', '100km', '1000km', 'LEO','GEO', 'Moon', r'10$^6$km',r'10$^7$km', 'Mercury', 'Venus', 'Earth','Mars', 'Jupiter', 'Satrun', 'Uranus', 'Neptune', 'Heliopause', 'Oorth', 'Outer Oorth', 'Alpha Centauri']

distLines = np.array([1e+1, 1e+2,   1e+3, 1e+4,   1e+5,    1e+6,     1e+7,       1e+8,        1e+9,        1e+10,       1.49597e+11,    1.49597e+12,  1.495979E+013,     0.94605284e+14,      0.94605284e+15,  0.94605284e+16, 0.94605284e+17 ])
distTexts =          ['10m','100m', '1km','10km', '100km', '1000km', r'10$^4$km',r'10$^5$km', r'10$^6$km', r'10$^7$km', '1AU',          '10AU',       '100AU',           '0.01ly',            '0.1ly',         '1ly',          '10ly'         ]


def plotShips():
    # --- setup
    times = 10**np.linspace( 0, 11.0, 100  )
    nodes = np.linspace( 0.0,  1.0, 100 ); # nodes = nodes**2; nodes[0] = NaN
    # big range
    tmin=times[0]; tmax=times[-1]
    dmin=1       ;  dmax=1e+17
    # --- main
    plt.figure(figsize=(15,10))
    plt.legend( loc=2, prop={'family': 'monospace'})
    for ship in ships:
	    print( ship[5], ship[1]*1e-3," km/s",ship[0]*1e-3, "kW/kg" )
	    st,vt = timeToDistace( times, ship[0], ship[1],     ship[2], ship[3], ship[4], nodes=nodes )
	    plt.plot( st, times, label=ship[5], color=ship[6], ls=ship[7] );
	    #text ( st[-1], times[-1], ship[5], color='k', horizontalalignment='left', verticalalignment='center')
    for time,txt in zip( timeLines, timeTexts ):
	    if ( time < tmax) and (time > tmin):
		    plt.axhline( time, color='k', alpha=0.5 )
		    plt.text   ( dmax, time, txt, color='k', horizontalalignment='right', verticalalignment='baseline' )
    for dist,txt in zip( distLines, distTexts ):
	    if ( dist < dmax) and (dist > dmin):
		    plt.axvline( dist, color='k', alpha=0.5 )
		    plt.text   ( dist, tmin, txt, rotation=90, color='k', horizontalalignment='left', rotation_mode='anchor', verticalalignment='baseline')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(dmin,dmax)
    plt.ylim(tmin,tmax)
    plt.ylabel(r"Flight Time [ s ]")
    plt.xlabel(r"Distance    [ m ]")
    plt.grid()
    #plt.savefig( "timeToDistance.png", bbox_inches='tight' )	
    plt.show()

#if __name__ == "__main__":

