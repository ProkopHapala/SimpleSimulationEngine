#from pylab import *

import numpy as np
import matplotlib.pyplot as plt

model_params={
'CD0'    : 0.008,  
'dCD'    : 0.5,  
'dCDS'   : 0.9,  
'dCL'    : 6.28,
'dCLS'   : 2.70,
'sStall' : 0.15,
'wStall' : 0.10,
}


# ==== Functions

def model_simple( sa, ca ):
	CD = np.abs(sa)
	CL = 2*sa*ca
	return CD,CL,CD*0


def model( sa, ca, params ):
	abs_sa = np.abs(sa)
	abs_ca = np.abs(ca)
	mask0  = abs_sa <   params['sStall']
	mask1  = abs_sa > ( params['sStall'] + params['wStall'] )
	mask01 = np.logical_not( np.logical_or( mask0 , mask1 ) )
	wS            = np.zeros(len(sa))
	wS[ mask1  ]  = 1
	a             = (abs_sa[mask01] - params['sStall'])/params['wStall']
	wS[ mask01 ]  = a*a*( 3 - 2*a )
	CD = params['CD0'] + ( (1-wS)*params['dCD']*abs_sa + wS*params['dCDS']    ) * abs_sa
	CL =                 ( (1-wS)*params['dCL']*ca     + wS*params['dCLS']*ca ) * sa 
	return CD,CL,wS

def multdata( F, sgn=1 ):
	n = len(F)
	F_=np.zeros( n*4 )
	F_[:n     ] = F
	F_[n:2*n  ] = F[::-1] *sgn
	F_[2*n:3*n] = F 
	F_[3*n:   ] = F[::-1] *sgn
	return F_

# ==== Main

# data
data = np.transpose(np.genfromtxt( 'NACA0008_Re 0.30_M0.00_N3.0.txt', skip_header=11 ))
CL = data[1];  CL    = multdata( CL, sgn=-1 ) 
CD = data[2];  CD    = multdata( CD, sgn= 1 )

# model
alfas = np.linspace( -np.pi,np.pi,len(CD) )
sa    = np.sin(alfas)
ca    = np.cos(alfas)
CD_,CL_,wS = model( sa, ca, model_params )
CD__,CL__,wS = model_simple( sa, ca )

#plot_1  angle dependence
plt.figure()
plt.plot( alfas, CD , 'r-', label='CD')
plt.plot( alfas, CL , 'b-', label='CL')
plt.plot( alfas, CD_, 'm-', label='CD_')
plt.plot( alfas, CL_, 'c-', label='CL_')
plt.plot( alfas, CD__, 'm--', label='CD__')
plt.plot( alfas, CL__, 'c--', label='CL__')
plt.plot( alfas, wS, 'k-', label='wS')
#axis('equal'); 
plt.grid(); 
plt.legend()

#plot_2  poloar_plot
plt.figure()
plt.plot( CD, CL  , 'k-' )
plt.plot( CD_, CL_, 'r-' )
plt.plot( CD__, CL__, 'r--' )
#plt.axis('equal'); 
plt.grid(); 
#plt.xlim(0.0,1.5)
	
plt.show()
