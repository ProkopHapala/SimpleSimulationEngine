#!/usr/bin/python

from pylab import *
import numpy as np
import os 


print "========= Load Library ============="

def makeclean():
	import os
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".so") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".o") ]
	[ os.remove(f) for f in os.listdir(".") if f.endswith(".pyc") ]

os.chdir('../');       print " >> WORKDIR: ", os.getcwd()
makeclean( )
sys.path.insert(0, "./")
import KosmoSuiteCpp as KS
os.chdir('examples');  print " >> WORKDIR: ", os.getcwd()


print "========= setup ============="

n = 10000;

ts   = linspace( 0, 10e-6, n );
dt   = ts[1] - ts[0] 
buff = zeros((n,7))

figure();

KS.FissionPulse_set_params ( generation_rate = 1/1e-7 );

v0s   = [ -1e+4, -2e+4, -5e+4, -10e+4 ]
clrs  = [  'b' , 'm',   'r',  'y'     ]


for i,v0 in enumerate(v0s):
	#KS.FissionPulse_set_initial( Nf0 = (8.0/0.239)*6.02214085774e+23, Nn0 = 1e+6, v0 = v0 );
	KS.FissionPulse_set_initial( Nf0 = (8.0/0.239)*6.02214085774e+23, Nn0 = 0, v0 = v0 );
	KS.FissionPulse_run_fixStep( buff, dt = dt );
	subplot(2,4,1); plot(  ts*1e+6, buff[ :, 0 ] , color=clrs[i] );  
	subplot(2,4,2); plot(  ts*1e+6, buff[ :, 1 ] , color=clrs[i] ); 
	subplot(2,4,3); plot(  ts*1e+6, buff[ :, 2 ] , color=clrs[i] );  
	subplot(2,4,4); plot(  ts*1e+6, buff[ :, 3 ] , color=clrs[i] ); 
	subplot(2,4,5); plot(  ts*1e+6, buff[ :, 4 ] , color=clrs[i] ); 
	subplot(2,4,6); plot(  ts*1e+6, buff[ :, 5 ] , color=clrs[i] ); 
	subplot(2,4,7); plot(  ts*1e+6, buff[ :, 6 ] , color=clrs[i] ); 
	subplot(2,4,8); plot(  ts*1e+6, buff[ :, 4 ], label = " v0= %2.2e m/s "%v0 ); 
	KS.FissionPulse_set_initial( Nf0 = 0, Nn0 = 0, v0 = v0 );
	KS.FissionPulse_run_fixStep( buff, dt = dt );
	subplot(2,4,1); plot(  ts*1e+6, buff[ :, 0 ] , color=clrs[i], ls='--' );  
	subplot(2,4,2); plot(  ts*1e+6, buff[ :, 1 ] , color=clrs[i], ls='--' ); 
	subplot(2,4,6); plot(  ts*1e+6, buff[ :, 5 ] , color=clrs[i], ls='--' ); 



subplot(2,4,1); grid(); xlabel(r'time[$\mu s$]');  ylabel(' R  [m]'     ); title( " Core Radius            "); ylim( 0.0, 0.2 ); 
subplot(2,4,2); grid(); xlabel(r'time[$\mu s$]');  ylabel(' v  [m/s]'   ); title( " Expansion Velocity     ");
subplot(2,4,3); grid(); xlabel(r'time[$\mu s$]');  ylabel(' Nf [1]'     ); title( " Fuel Nuclei population "); 
subplot(2,4,4); grid(); xlabel(r'time[$\mu s$]');  ylabel(' Nn [1]'     ); title( " Neutron population     "); yscale('log')
subplot(2,4,5); grid(); xlabel(r'time[$\mu s$]');  ylabel(' Q  [J]'     ); title( " Total Energy           "); yscale('log')
subplot(2,4,6); grid(); xlabel(r'time[$\mu s$]');  ylabel(' F  [N]'     ); title( " Plasma pressure Force  ");
subplot(2,4,7); grid(); xlabel(r'time[$\mu s$]');  ylabel(' P_inter [1]');
subplot(2,4,8); legend();

#xs = linspace(1,10,100)
#plot( xs, sin(xs) ); grid()


show()
 
