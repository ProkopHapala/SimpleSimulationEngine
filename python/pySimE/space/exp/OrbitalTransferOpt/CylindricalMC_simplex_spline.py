#!/usr/bin/env python

from pylab import *
from Poly6th_numeric import *

from Simplex2_prokop import Simplex
from CubicClampedSpline import *

T= 3.0

nnodes  = 4
perNode = 10

ts = array( range( (nnodes+1)*perNode )  )*T

print " len(ts) : ", len(ts)

def CircularP0V0( R0, T, phi0, vr ):
    cdO0  = sqrt(1/R0**3)
    #return array([  R0 ,  cdO0*phi0*T  ]), array( [ vr, cdO0 ] ) * T
    return array([  R0 ,  2*pi*phi0  ]), array( [ vr, cdO0 ] )
    
P0,V0 = CircularP0V0( 1.0, T, 0.0, 0 )
P1,V1 = CircularP0V0( 0.2, T, 2  , 0 )

def initGen(P0,P1,n):
    Gen= [0.0]*(2*n)
    for i in range(n):
        dt=float(i+1)/(n+2); mt= 1.0-dt
        Gen[i  ]= mt*P0[0] + dt*P1[0] 
        Gen[i+n]= mt*P0[1] + dt*P1[1]
    return Gen

def evalGen ( Gen ):
    n = (len(Gen)/2)
    ti = [0]*(n+2); ri = [0]*(n+2); oi = [0]*(n+2);
    ti[0  ]=0; ri[0]  =P0[0]; oi[0]  =P0[1];
    ti[n+1]=T; ri[n+1]=P1[0]; oi[n+1]=P1[1];
    for i in range(0,n):
        ti[i+1]=(i+1)*T/(n+1)
        ri[i+1]=Gen[i  ]
        oi[i+1]=Gen[i+n]
    #print len(ti),len(ri)
    SR, arn = Spline4_clamped( ti, ri, V0[0], V1[0])
    SO, aon = Spline4_clamped( ti, oi, V0[1], V1[1])
    #print " len(SR) ",len(SR)," len(SO) ", len(SO)
    #CRs = poly6Coefs( P0[0], V0[0],  Gen[0],   P1[0], V1[0], Gen[1] )
    #COs = poly6Coefs( P0[1], V0[1],  Gen[2],   P1[1], V1[1], Gen[3] ) 
    return evalPolarForces( T, SR, SO, perNode)

maxThrust = 2.0
def fitnesFunc( Fs ):
    fsum  = 0
    tsum = 0
    #print "len(Fs) ", len(Fs[4])," len(ts) ", len(ts)
    for i in range(len(Fs[4])-1):
        dt=(ts[i+1]-ts[i])
        #df=0.5*(Fs[4][i+1]+Fs[4][i])
        df=0.5*(Fs[4][i+1]**2+Fs[4][i]**2)
        fsum+=df*dt
        tsum+=dt
        #df_over = df-maxThrust
        #if(df_over>0):
        #    fsum+= (df_over**2) * dt   # penalty for overloading engine
    return -fsum
    #return -T*  sqrt((Fs[4]**2).sum()) /len(ts)
    
def evalFitness( Gen ):
    global Os,Rs,Fs
    Os,Rs,Fs   =  evalGen ( Gen )
    fitness =  fitnesFunc(Fs)
    return -fitness

def plotTrj( Os,Rs,Fs, i ):
    #subplot(2,5,1+5*i); plot( Os[0], Rs[0], '.-' );  grid()
    subplot(2,5,1+5*i, polar=True); plot( Os[0], Rs[0], '.-');  
    subplot(2,5,2+5*i);  plot( ts, Rs[0],'r-' ); plot( ts, Os[0], 'b-' );grid()
    subplot(2,5,3+5*i);  plot( ts, Rs[1],'r-' ); plot( ts, Os[1], 'b-' );grid() 
    subplot(2,5,4+5*i); 
    plot( ts, Rs[2],'r:' );  plot( ts, Os[2], 'b:' );   
    plot( ts, Fs[1],'r-' );  plot( ts, Fs[0], 'b.-' );   
    plot( ts, Fs[2],'r--' );   # G
    plot( ts, Fs[3],'r.-' );   # FTR
    plot( ts, Fs[4],'k.-' );   # FT 
    grid()

# ================ MAIN PROGRAM BODY =========================
figure(num=None, figsize=(20, 10))

GenHistory = []

#Gen = array([ A0es[0], A1es[0], A0es[1], A1es[1] ])
#Gen = MCrun( Gen, 10.0,  0.001, 100, 10  )
Gen = initGen(P0,P1, nnodes)
print " Initial Gen : ", Gen

Os,Rs,Fs = evalGen(Gen)

print "shape(ts)", shape(ts)
print "shape(ts)", shape(Fs[4])
#fitness = fitnesFunc(Fs)
plotTrj( Os,Rs,Fs, 0 )

#Simp = Simplex(evalFitness, Gen, [1, 1, 1, 1])

steps = ones(nnodes*2)*0.1
Simp = Simplex(evalFitness, Gen, steps )
#values, err, niter = SimplexOpt.minimize()
old_low = 10000000000
for i in range(0, 1000):
    converged, err,low,hi = Simp.simplexStep( 0.01 )
    if converged:
        print " converged in ",i," steps " 
        break;
    if(low < old_low):
        old_low = low 
        subplot(2,5,5); plot( ts, Fs[4], '-', lw=0.25 );  grid()
        GenHistory.append(list(Simp.simplex[Simp.lowest]))
        print " new_low :  ", low, " iter: ", i, " err ", err
print Simp.simplex[Simp.lowest]
subplot(2,5,5); plot( ts, Fs[4], 'k.-', lw=2 );  grid()

Os,Rs,Fs = evalGen( Gen)
plotTrj( Os,Rs,Fs, 1 )


if len(GenHistory)>2:
    GenHistory = transpose(array(GenHistory ))
    subplot(2,5,10); 
    plot( GenHistory[0], '.-', label="ar0" ); 
    plot( GenHistory[1], '.-', label="ar1" ); 
    plot( GenHistory[2], '.-', label="ao0" ); 
    plot( GenHistory[3], '.-', label="ao1" ); 
    legend( bbox_to_anchor=(0.5, 1.00, 1., 0.000) )


show()