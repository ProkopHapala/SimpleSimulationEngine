import matplotlib.pyplot as plt
import numpy as np

# =================== Functions

def tsielkovsky_deltav( vexh, mass_ratio ):
	return vexh * np.log( mass_ratio )
	
# =================== deltaV-vs-trust

times = np.array([1.0,  60.0, 3600.0, 86400.0, 604800,  2592000, 31556926, 315569260,  3155692600, 31556926000 ])
timeText = ['sec','min','hour', 'day',   'week',  'month', 'year',   '10years',  '100years', '1000years']

dists= np.array([1e+1, 1e+2, 1e+3, 1e+4, 1e+5, 1e+6, 6371e+3, 42164e+3, 384400e+3, 1e+9, 1e+10, 5.790918E+010, 1.082089E+011, 1.495979E+011, 2.279366E+011, 7.784120E+011, 1.426725E+012, 2.870972E+012,  4.498253E+012, 1.40621998e+13, 2.99195741e+14, 7.47989354e+15, 4.13425091e+16])
distText      = ['10m','100m', '1km','10km', '100km', '1000km', 'LEO','GEO', 'Moon', r'10$^6$km',r'10$^7$km', 'Mercury', 'Venus', 'Earth','Mars', 'Jupiter', 'Satrun', 'Uranus', 'Neptune', 'Heliopause', 'Oorth', 'Outer Oorth', 'Alpha Centauri']

powerDens = np.array([ 1.0,    10,      100,       1e+3,       1e+4,      1e+5,       1e+6,       1e+7,    1e+8       ])
powerText =          ['1W/kg','10W/kg','100W/kg', '1kW/kg',  '10kW/kg',  '100kW/kg', '1MW/kg',  '10MW/kg', '100MW/kg' ]


'''
#  name                delta-v[km/s] accel[m/s^2]
ships=[
['Boeing1981',			51082.56,	0.01600  ],
['NERVA',				13750.33,	9.80628  ],
#['Bimodal NTR',			11458.58,	8.49257  ],
#['Liberty Bell',		3917.25,	15.45000 ],
['Discovery2',			225258.04,	0.02039  ],
#['ExactingStarfighter',	6591673.73,	150.00000],
#['LibertyShip',			15697.44,	23.36250 ],
['HELIOS',				20968.33,	144.26471],
['HydeFussion',			764401.46,	0.02058  ],
['HOPE_FFRE',			138336.01,	0.00015  ],
['HOPE_MPD',			22367.09,	0.00006  ],
#['Stuhlinger_Ion',		165495.96,	0.00041  ],
['HOPE_Zpinch',			55214.14,	0.06906  ],
['ICAN-2',				94708.18,	0.52174  ],
['LCOTV',				8093.51,	0.00051  ],
['MarsUmberla',			58685.9,	0.00069  ],
['MiniOrionDRM-2',		99790.9,	2.56493  ],
#['IBS Agamemnon',		280052.45,	0.35714  ],
#['Orion',				94494.94,	87.91209 ],
['VISTA',				201401.55,	0.13079  ]
]
'''

#  name                delta-v[km/s] accel[m/s^2]
ships=[
    ['NERVA',				13750.33,	9.80628  ],
    ['LCOTV_Solar_Ion',				8093.51,	0.00051  ],

    ['Discovery2' ,			225258.04,	0.02039  ],
    ['HOPE_FFRE'  ,			138336.01,	0.00015  ],
    ['HOPE_Zpinch',			55214.14,	0.06906  ],
    #['MiniMagOrion'  ,		99790.9,	2.56493      ],

    ['MiniMagOrion'  ,		tsielkovsky_deltav( 93e+3, 788686.0/157723.0  ),	642e+3/( (788686+157723) * 0.5 )      ],

    ['Orion_10kt',      tsielkovsky_deltav( 120e+3, 10e+6/4.348e+6  ),      80e+6/( ( 10e+6 + 4.348e+6) * 0.5 ) ],

    ['Orion_100m',      tsielkovsky_deltav( 10.0e+6, 400e+3/100e+3  ), 9.81 ],
    ['Daedalus_1st',   tsielkovsky_deltav( 10.6e+6, 46e+6/1690e+3  ),     7.540e+6/(46e+6*0.5)       ],
]


def main():
    vMin, vMax = 1e+0, 1e+8
    aMin, aMax = 1e-5, 1e+3
    As = np.linspace(aMin, aMax, 2)
    fig, ax = plt.subplots(figsize=(12, 8))
    for dist, text in zip(dists, distText):
	    # compute the line
	    v = np.sqrt(2*dist * As)
	    ax.plot(v, As, 'b-', alpha=0.5)
	    # sort out where the label should be
	    txt_y = aMin
	    txt_x = v[0]
	    # clip to the edges
	    if (txt_x < vMin):
		    txt_x = vMin
		    txt_y = vMin**2 / (2*dist)
	    ax.text(txt_x, txt_y, text,  rotation=55, color='b', rotation_mode='anchor', horizontalalignment='left', verticalalignment='bottom')
    for time, txt in zip(times, timeText):
	    # compute the line
	    v = As * time
	    ax.plot(v, As, 'r-', alpha=0.5) 
	    # sort out where the label should be
	    txt_x = v[-1]
	    txt_y = aMax
	    # clip to the edges
	    if(txt_x > vMax):
		    txt_x = vMax
		    txt_y = vMax / time
	    ax.text(txt_x, txt_y, txt, rotation=35, color='r', horizontalalignment='right', rotation_mode='anchor', verticalalignment='baseline')
    for power, txt in zip(powerDens, powerText):
	    # compute the line
	    v = power / As
	    ax.plot(v, As, 'g-', alpha=0.5) 
	    # sort out where the label should be
	    txt_x = v[-1]
	    txt_y = aMax
	    # clip to the edges
	    if(txt_x < vMin):
		    txt_x = vMin
		    txt_y = power / vMin
	    ax.text(txt_x, txt_y, txt, rotation=-40, color='g', horizontalalignment='left', rotation_mode='anchor', verticalalignment='top')
    for ship in ships:
	    ax.plot( ship[1], ship[2], 'ok' )
	    ax.text( ship[1], ship[2], ship[0], color='k', rotation=0, rotation_mode='anchor' )
    ax.set_ylabel(r"acceleration [m/s$^2$]")
    ax.set_xlabel(r"delta-v [m/s]")
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid()
    ax.set_ylim(aMin, aMax)
    ax.set_xlim(vMin, vMax)
    #plt.savefig( 'deltaV_accel_limits.png', bbox_inches='tight' );
    plt.show()
    
if __name__ == "__main__":
    main()
