from pylab import *

def model( sa, ca, params ):
	abs_sa = abs(sa)
	abs_ca = abs(ca)
	mask0  = abs_sa <   params['sStall']
	mask1  = abs_sa > ( params['sStall'] + params['wStall'] )
	mask01 = logical_not( logical_or( mask0 , mask1 ) )
	wS            = zeros(len(sa))
	wS[ mask1  ]  = 1
	a             = (abs_sa[mask01] - params['sStall'])/params['wStall']
	wS[ mask01 ]  = a*a*( 3 - 2*a )
	CD = params['CD0'] + ( (1-wS)*params['dCD']*abs_sa + wS*params['dCDS']        ) * abs_sa
	CL =                 ( (1-wS)*params['dCL']        + wS*params['dCLS']*abs_ca ) * sa 
	return CD,CL,wS


alfas = linspace( -pi,pi,1000 )

'''
paramsRudder={
'CD0'    : 0.008,  
'dCD'    : 1.5,  
'dCDS'   : 0.9,  
'dCL'    : 6.28,
'dCLS'   : 2.70,
'sStall' : 0.16,
'wStall' : 0.08,
}

paramsSail={
'CD0'    : 0.1,  
'dCD'    : -1.0,  
'dCDS'   : 0.8,  
'dCL'    : 3.50,
'dCLS'   : 2.20,
'sStall' : 0.20,
'wStall' : 0.40,
}

paramsKeel={
'CD0'    : 0.04,  
'dCD'    : 1.5,  
'dCDS'   : 0.9,  
'dCL'    : 3.00,
'dCLS'   : 2.00,
'sStall' : 0.20,
'wStall' : 0.40,
}
'''


paramsRudder={
'CD0'    : 0.008,  
'dCD'    : 0.9,  
'dCDS'   : 0.9,  
'dCL'    : 2.70,
'dCLS'   : 2.70,
'sStall' : 0.16,
'wStall' : 0.08,
}

paramsSail={
'CD0'    : 0.1,  
'dCD'    : 0.8,  
'dCDS'   : 0.8,  
'dCL'    : 2.20,
'dCLS'   : 2.20,
'sStall' : 0.20,
'wStall' : 0.40,
}

paramsKeel={
'CD0'    : 0.04,  
'dCD'    : 0.9,  
'dCDS'   : 0.9,  
'dCL'    : 2.00,
'dCLS'   : 2.00,
'sStall' : 0.20,
'wStall' : 0.40,
}


sa  = sin(alfas)
ca  = cos(alfas)


CD_R,CL_R,wS_R = model( sa, ca, paramsRudder );    print " rudder max(L/D) : ", max(CL_R/CD_R)
CD_S,CL_S,wS_S = model( sa, ca, paramsSail   );    print " Sail   max(L/D) : ", max(CL_S/CD_S)
CD_K,CL_K,wS_K = model( sa, ca, paramsKeel   );    print " Keel   max(L/D) : ", max(CL_K/CD_K)

'''
figure()
plot( alfas, CD_R , 'r-', label='CD')
plot( alfas, CL_R , 'b-', label='CL')
plot( alfas, wS_R, 'k-', label='wS')

plot( alfas, CD_S , 'r-', label='CD')
plot( alfas, CL_S , 'b-', label='CL')
plot( alfas, wS_S, 'k-', label='wS')

plot( alfas, CD_K , 'r-', label='CD')
plot( alfas, CL_K , 'b-', label='CL')
plot( alfas, wS_K, 'k-', label='wS')

axis('equal'); grid(); legend()
'''




figure()
plot( alfas, CL_R, 'r-', label='Rudder' )
plot( alfas, CL_S, 'g-', label='Sail')
plot( alfas, CL_K, 'b-', label='Keel')
grid(); legend()


figure()
plot( alfas, CD_R, 'r-', label='Rudder' )
plot( alfas, CD_S, 'g-', label='Sail')
plot( alfas, CD_K, 'b-', label='Keel')
grid(); legend()

figure()
plot( CD_R, CL_R, 'r-', label='Rudder' )
plot( CD_S, CL_S, 'g-', label='Sail')
plot( CD_K, CL_K, 'b-', label='Keel')
grid(); legend()
	
show()
