# -*- coding: utf-8 -*-
"""
Created on Fri May 03 20:54:06 2013

@author: asiJa
"""

from pylab import *

N  = 127
dt = 1.0/N 

ts = arange(0,1.000001,dt)

def evalCos(nf,ts):
	f = pi*nf;
	cs = cos(f*ts)
	ss = sin(f*ts)
	Y   =      cs
	dY  =   -f*ss
	ddY = -f*f*cs
	return Y,dY,ddY

def evalDiCos(nf1,nf2, ts):
	f1  = pi*nf1; 	      f2 = pi*nf2 
	cs1 = cos(f1*ts); 	cs2 = cos(f2*ts)
	ss1 = sin(f1*ts); 	ss2 = sin(f2*ts)
	Y   =        cs1 -     cs2
	dY  =    -f1*ss1 +     f2*ss2
	ddY = -f1*f1*cs1 +     f2*f2*cs2
	return Y,dY,ddY

def evalDiCosTransf(V,ts):
	n = len(ts)
	Y=zeros(n);dY=zeros(n);ddY=zeros(n);
	for i in range(len(V)):
		C = V[i]/((i+1)**2)
		f1  = pi*i; 	      f2 = pi*(i+2) 
		cs1 = cos(f1*ts); 	cs2 = cos(f2*ts)
		ss1 = sin(f1*ts); 	ss2 = sin(f2*ts)
		Y   +=C*(        cs1 -    cs2       )
		dY  +=C*(    -f1*ss1 +    f2*ss2    )
		ddY +=C*(  -f1*f1*cs1 +   f2*f2*cs2 )
	return Y,dY,ddY
	
#n1=1;n2=3;
#n1=2;n2=4;

'''
Y1,dY1,ddY1 = evalCos(n1, ts)
Y2,dY2,ddY2 = evalCos(n2, ts)
Y,dY,ddY    = evalDiCos(n1,n2, ts)

subplot(1,3,1); plot(ts,  Y,"-k"); plot(ts,  Y1,"-b"); plot(ts,  -Y2,"-r"); 
subplot(1,3,2); plot(ts, dY,"-k"); plot(ts, dY1,"-b"); plot(ts, -dY2,"-r"); 
subplot(1,3,3); plot(ts,ddY,"-k"); plot(ts,ddY1,"-b"); plot(ts,-ddY2,"-r"); 
'''

'''
figure(num=None, figsize=(14, 8))
for i in range(3):
	Y,dY,ddY = evalDiCos(i,i+2, ts)
	subplot(1,3,1); plot(ts,  Y);
	subplot(1,3,2); plot(ts, dY);
	subplot(1,3,3); plot(ts,ddY);
'''

'''
figure(num=None, figsize=(14, 8))
for i in range(0,16,4):
	Y,dY,ddY = evalDiCos(i,i+2, ts)
	subplot(1,3,1); plot(ts,  Y);
	subplot(1,3,2); plot(ts, dY);
	subplot(1,3,3); plot(ts,ddY);
	Y,dY,ddY = evalDiCos(i+1,i+3, ts)
	subplot(1,3,1); plot(ts,  Y);
	subplot(1,3,2); plot(ts, dY);
	subplot(1,3,3); plot(ts,ddY);
'''

V = [0.01,0.04,-0.03,-0.2]

Y,dY,ddY  = evalDiCosTransf(V,ts)

subplot(1,3,1); plot(ts,  Y);
subplot(1,3,2); plot(ts, dY);
subplot(1,3,3); plot(ts,ddY);

show()