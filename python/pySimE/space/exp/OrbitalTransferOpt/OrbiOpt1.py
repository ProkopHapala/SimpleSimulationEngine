from pylab import *
from Poly6th_numeric import *


# ==== boundary condition on position and force
P0 = array([  0   ,  0  , 0]);   V0 = array([ 0.5, 1.0, 0])  
P1 = array([  1.0 ,1.0  , 0]);   V1 = array([ 0.0, 0.5, 0])
 
# ==== Estimate Acceleleration
#A0es,A1es = ForceEstimate_Linear ( P0, V0, P1, V1 )
A0es,A1es = ForceEstimate_Cos       ( P0, V0, P1, V1 )
print " A0es :", A0es, " A1es :", A1es

A0 = A0es
A1 = A1es

print " A0   :", A0  , " A1   :", A1

As = poly6Coefs( P0[0], V0[0],  A0[0],   P1[0], V1[0], A1[0] )
print " As :", As
print " P0 V0 A0: ", P0[0], V0[0], A0[0], " Poly(t=0) : ",   evalPoly6( 0.0, As)
print " P0 V0 A0: ", P1[0], V1[0], A1[0], " Poly(t=1) : ",   evalPoly6( 1.0, As)


# Evaluate Arrays
ts = arange(0,1.0001,0.01)
xs = zeros([3,len(ts)])
vs = zeros([3,len(ts)])
fs = zeros([3,len(ts)])
for i in range(3):
    As = poly6Coefs( P0[i], V0[i],  A0[i],   P1[i], V1[i], A1[i] )
    evalPoly6_array( ts, As,   xs[i],vs[i],fs[i] )


# ======== JUST POTTING ==================

colors = [ 'r','g','b' ]
subplot(2,3,1); title(" positions ")
for i in range(3):
    plot(ts, xs[i], color = colors[i])
    
subplot(2,3,2);  title(" velocity ")
for i in range(3):
    plot(ts, vs[i], color = colors[i])

subplot(2,3,3);  title(" acceleration ")
for i in range(3):
    plot(ts, fs[i], color = colors[i])  

subplot(2,3,4);   plot(xs[0], xs[1], '-k')
subplot(2,3,5);   plot(vs[0], vs[1], '-k')  
subplot(2,3,6);   plot(fs[0], fs[1], '-k') 

show()