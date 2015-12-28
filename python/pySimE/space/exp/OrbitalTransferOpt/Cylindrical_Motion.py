#!/usr/bin/env python

from pylab import *
from Poly6th_numeric import *


def CircularP0V0( R0,T, t ):
    cdO0  = sqrt(1/R0**3)
    return array([  R0 ,  cdO0*t  ]), array( [ 0.0, cdO0 ] )*T

T= 2.0
P0,V0 = CircularP0V0( 1.3, T, 0.0 )
P1,V1 = CircularP0V0( 0.5, T, T/2   )

# ==== boundary condition on position and force
#P0 = array([  cR0 ,  0  ]);   V0 = array( [ 0.0, cdO0 ] )*T  
#P1 = array([  cR0 ,  1  ]);   V1 = array( [ 0.0, cdO0 ] )*T
 
# ==== Estimate Acceleleration
#A0es,A1es = ForceEstimate_Linear ( P0, V0, P1, V1 )
A0es,A1es = ForceEstimate_Cos       ( P0, V0, P1, V1 )
print " A0es :", A0es, " A1es :", A1es

A0 = A0es
A1 = A1es

# Evaluate Arrays
ts = arange(0,1.0001,0.05)
Rs = poly6Coefs( P0[0], V0[0],  A0[0],   P1[0], V1[0], A1[0] )
Os = poly6Coefs( P0[1], V0[1],  A0[1],   P1[1], V1[1], A1[1] )
Os,Rs,Fs = evalPolarForces(ts, T,  Rs, Os)

#print Rs[0]
#print Os[0]


print " R    :   ", Rs  [0]
print " dO   :   ", Os  [1]
print " ddR  :   ", Rs  [2]
print " F_R  :   ", Fs  [1]

figure(num=None, figsize=(16, 4))

subplot(1,4,1 ); plot( Os[0], Rs[0], '.-' );  axis('equal'); grid()
#subplot(4,3,2 ); plot( Os[1], Rs[1], '.-' );  axis('equal'); grid()
#subplot(4,3,3 ); plot( Fs[0], Fs[1], '.-' );  axis('equal'); grid()

subplot(1,4,2); plot( ts, Rs[0],'r-' ); plot( ts, Os[0], 'b-' );grid()
subplot(1,4,3); plot( ts, Rs[1],'r-' ); plot( ts, Os[1], 'b-' );grid() 
subplot(1,4,4); 
plot( ts, Rs[2],'r:' );  plot( ts, Os[2], 'b:' );   
plot( ts, Fs[1],'r-' );  plot( ts, Fs[0], 'b.-' );   
plot( ts, Fs[2],'r--' );   # G
plot( ts, Fs[3],'r.-' );   # FTR
plot( ts, Fs[4],'k.-' );   # FT 
grid()




'''
coss = cos(Os[0])
sins = sin(Os[0])
rhat = [ coss, sins ]
ohat = [-sins, coss ]

xs  = Rs[0]*rhat[0];   ys = Rs[0]*rhat[1]
vxs = Rs[1]*rhat[0] +  Os[1]*ohat[0];
vys = Rs[1]*rhat[1] +  Os[1]*ohat[1];
#axs = Fs[]
subplot(3,3,1); plot(xs, ys );    axis('equal');
subplot(3,3,2); plot(vxs, vys );  axis('equal');
'''


# ======== JUST POTTING ==================


'''
subplot(3,3,1); title(" positions ")
plot(ts, xs[i], color = colors[i])
    
subplot(3,3,2);  title(" velocity ")
plot(ts, vs[i], color = colors[i])

subplot(3,3,3);  title(" acceleration ")
for i in range(3):
plot(ts, fs[i], color = colors[i])  

subplot(2,3,4);   plot(xs[0], xs[1], '-k')
subplot(2,3,5);   plot(vs[0], vs[1], '-k')  
subplot(2,3,6);   plot(fs[0], fs[1], '-k') 
'''

show()