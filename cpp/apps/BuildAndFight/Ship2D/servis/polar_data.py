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


data = transpose(genfromtxt( 'NACA0008_Re 0.30_M0.00_N3.0.txt', skip_header=11 ))

CL = data[1]
CD = data[2]

def multdata( F, sgn=1 ):
	n = len(F)
	F_=zeros( n*4 )
	F_[:n     ] = F*sgn
	F_[n:2*n  ] = F[::-1]*sgn 
	F_[2*n:3*n] = F 
	F_[3*n:   ] = F[::-1]
	return F_



CD = multdata( CD, sgn=1 )
CL = multdata( CL, sgn=-1 )
alfas = linspace( -pi,pi,len(CD) )
sa  = sin(alfas)
ca  = cos(alfas)

params={
'CD0'    : 0.008,  
'dCD'    : 0.5,  
'dCDS'   : 0.9,  
'dCL'    : 6.28,
'dCLS'   : 2.70,
'sStall' : 0.16,
'wStall' : 0.08,
}

CD_,CL_,wS = model( sa, ca, params )

figure()
plot( alfas, CD , 'r-', label='CD')
plot( alfas, CL , 'b-', label='CL')
plot( alfas, CD_, 'r--', label='CD_')
plot( alfas, CL_, 'b--', label='CL_')
plot( alfas, wS, 'k-', label='wS')

#axis('equal'); 
grid(); 
legend()

figure()
plot( CD, CL  , 'k-' )
plot( CD_, CL_, 'r-' )
#axis('equal'); 
grid(); 
#xlim(0.0,1.5)
	
show()
