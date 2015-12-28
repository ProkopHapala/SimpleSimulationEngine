# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 08:35:21 2013

@author: asiJa
"""
from pylab import *
from scipy.integrate import odeint
from scipy.optimize import minimize_scalar


specific_strength={ 'S-Glass':1988e+3, 'Basalt':1790e+3, 'Vectran':2071e+3, 'Carbon':2071e+3, 'Kevlar':2514e+3, 'Dyneema':3711e+3, 'Zylon':3766e+3, 'PBOmax':5843e+3 , 'CNT': 46268e+3 }

def skyHookTaper( rTouch=6378.1e+3, lrope=1000e+3, vTouch=1400, s=3766e+3, M=5.9736e+24, G = 6.67384e-11 ):
	rSat = rTouch+lrope
	vSat     = sqrt ( G*M / rSat )
	r0       = rSat - lrope / ( 1.0 - vTouch/vSat )
	l0       = rSat - r0  
	#print " lrope,r0-rSat ", lrope/1000,rSat-r0/1000 
	omega    = vSat/rSat
	xs       = arange(0.0,lrope+0.0001, lrope/1000.0) 
	'''
	phis       = arange(0.0,2*pi+0.0001, pi/50.0)
	plot(rSat*cos(phis),rSat*sin(phis),'--', label='orbit')
	plot(r0*cos(phis),r0*sin(phis),':', label='r_v0')
	plot(rTouch*cos(phis),rTouch*sin(phis), label='r_ground' )
	def Epitrochoid(phi, r, h, x ):
		xs = r * cos(phi) - x*cos(phi*r/h)
		ys = r * sin(phi) - x*sin(phi*r/h)
		return xs,ys
	xs,ys =  Epitrochoid(phis, rSat, rSat-r0, lrope ); plot(xs,ys, '-k',lw=2,label='trajectory')
	xs,ys =  Epitrochoid(phis, rSat, rSat-r0, rSat-r0 ); plot(xs,ys, ':', label='v0 cycloide')
	axis('equal')	
	figure()
	aCs   =  omega**2  *  rSat*( (lrope-xs)*rSat/l0**2 - 1)  # accordign to x-acceleration oc cycloide at point phi=0 (second derivative of cartesian)
	aGs = M*G/(rTouch+xs)**2
	plot(xs,aCs,'.-r' )
	plot(xs,aGs,'.-g' )
	plot(xs,aCs+aGs,'.-b' )
	plot(xs,zeros(len(xs)) )
	'''
	def func( m, x, s=s, rTouch=rTouch, lrope=lrope,  C1=M*G, C2=omega**2*rSat, C3=rSat/l0**2 ):
		#print lrope
		return (m/s) *( C1/(rTouch + x)**2 + C2*( (lrope-x)*C3 -1 ) )
	ms = odeint( func, 1.0, xs, args=(s, rTouch, lrope,  M*G, omega**2*rSat, rSat/l0**2, ) )
	#figure(); 	plot(xs,ms )
	return xs , ms 

#skyHookTaper( s=specific_strength['Zylon'], lrope=100, vTouch=0 )


material        = 'Zylon'
strength_factor = 0.8
vTouch          = 1000
sstrength       = specific_strength[material]*strength_factor
titlestr=(" v0 = %1.1f km/s " %(vTouch/1000.0)) + material+' * '+str(strength_factor) ;
savestr =("v0=%1.1fkmps_" %(vTouch/1000.0)) + material+'x'+str(strength_factor) ;


mtots  = []
tapers = []
lropes = arange(0.1e+3,4000e+3,100e+3)
for i in range(len(lropes)):
	lrope= lropes[i]
	xs,ms = skyHookTaper( s=sstrength, lrope=lrope, vTouch=vTouch )
	m = trapz(transpose(ms), xs)
	mtots.append( m )
	tapers.append( ms[-1][0] )
	print ' lrope= ',lrope/1000.0,'[ km ] m=',m,' [kg] Taper= ' ,ms[-1][0]

subplot(2,1,1); plot( lropes/1000.0, array(mtots)/sstrength, '.-' );   xlabel('length [km]'); ylabel('mass'); 
title( titlestr);
subplot(2,1,2); plot( lropes/1000.0, tapers, '.-' );  xlabel('length [km]'); ylabel('taper');
savefig(savestr+'.png');


'''
material = 'Zylon'
strength_factor = 0.25

vTouchs = arange(0,3000,500)
xmins = []
ymins = []
for i in range(len(vTouchs)):
	mls = [];
	lropes = arange(100e+3,5000e+3,100e+3)
	def f(lrope):
		xs,ms = skyHookTaper( s=specific_strength[material]*strength_factor,lrope=lrope, vTouch=vTouchs[i] )
		m = trapz(transpose(ms), xs)
		#return ms[-1]
		return m
	for lrope in lropes:
		mls.append( f(lrope) ) 
	resmin = minimize_scalar(f, bounds=(250e+3,5000e+3), method='Bounded')
	xmins.append(   resmin.x   )
	ymins.append( f(resmin.x)  )
	print  ' vTouch= ',vTouchs[i]/1000.0, ' [km/s]  l_opt= ',xmins[-1]/1000.0, ' [km] Taper= ',ymins[-1]
	plot(lropes/1000,mls,'.-', label=('%1.1f km/s' %(vTouchs[i]/1000.0) ) ); 

plot(array(xmins)/1000,ymins,'.-', );
xlabel(' rope length [km]'); ylabel(' Taper [1]');   yscale('log'); legend()
title( " Sapce Rotovator of "+material+' * '+str(strength_factor) );
'''

show()