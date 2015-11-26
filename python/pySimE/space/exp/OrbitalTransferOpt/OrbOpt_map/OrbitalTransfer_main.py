#!/usr/bin/env python

from pylab import *

from Simplex_optimization import Simplex
from ClampedCubicSpline import *




nnodes   = 8
nsplines = nnodes + 1 
perNode  = 10
nsamp    = (nsplines)*perNode 
    
Gen= [0.0]*2*nnodes

def R2omega(R):
    return sqrt(1.0/R**3)


T  = 3.0
R0 = 1.0
R1 = 0.2
v0=R2omega(R0)
v1=R2omega(R1)

P0=array( [R0,0] );    V0=array( [1.5,v0] ); 

vo = v0
ph = 0.0 

dt = T/nsplines;

for i in range(nnodes):
    w=(i+1)/float(nsplines); mw= 1.0-w
    ri       = mw*R0 + w*R1
    Gen[i]   = ri
    vi       = R2omega(ri)
    ph      += 0.5*(vo+vi)*dt
    Gen[i+nnodes] = ph
    vo       = vi

vi  = R2omega(R1)
ph += 0.5*(vo+vi)*dt
P1=array( [R1,ph] );    V1=array( [-2.0,v1] ); 

print " P0 ", P0
print " P1 ", P1
print " V0 ", V0
print " V1 ", V1

ts = array( range(nsamp) )*T/float(nsamp)

def gen2controPoints( Gen ):
    n = (len(Gen)/2)
    ti = [0]*(n+2); ri = [0]*(n+2); oi = [0]*(n+2);
    ti[0  ]=0; ri[0]  =P0[0]; oi[0]  =P0[1];
    ti[n+1]=T; ri[n+1]=P1[0]; oi[n+1]=P1[1];
    for i in range(0,n):
        ti[i+1]=(i+1)*T/(n+1)
        ri[i+1]=Gen[i  ]
        oi[i+1]=Gen[i+n]
    return ti,ri,oi


def evalGen ( Gen ):
    #print len(ti),len(ri)
    ti,ri,oi =  gen2controPoints( Gen )
    n = len(ti)
    dt0=ti[1  ]-ti[0  ]
    dtn=ti[n-1]-ti[n-2]
    #SR, arn = Spline4_clamped( ti, ri, V0[0]/dt0, V1[0]/dtn)
    #SO, aon = Spline4_clamped( ti, oi, V0[1]/dt0, V1[1]/dtn)
    #print " len ti, ri,oi",len(ti),len(ri),len(oi)
    SR, arn = Spline4_clamped( ti, ri, V0[0], V1[0])
    SO, aon = Spline4_clamped( ti, oi, V0[1], V1[1])
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
    return -sqrt(fsum/tsum)
    #return -T*  sqrt((Fs[4]**2).sum()) /len(ts)
    
def evalFitness( Gen ):
    global Os,Rs,Fs
    Os,Rs,Fs   =  evalGen ( Gen )
    fitness =  fitnesFunc(Fs)
    return -fitness

def plotTrj( Os,Rs,Fs, i ):
    #subplot(2,5,1+5*i); plot( Os[0], Rs[0], '.-' );  grid()
    ti,ri,oi= gen2controPoints(Gen)
    subplot(2,5,1+5*i, polar=True); plot( Os[0], Rs[0], '.-k');  plot( oi, ri, 'ok');  
    subplot(2,5,2+5*i); 
    print len(ti),len(ri),len(oi)
    plot( ts, Rs[0],'r-' ); plot( ts, Os[0], 'b-' ); 
    plot( ti, ri,'ro' ); plot( ti, oi, 'bo' ); 
    grid()
    subplot(2,5,3+5*i);  plot( ts, Rs[1],'r-' ); plot( ts, Os[1], 'b-' ); grid() 
    subplot(2,5,4+5*i); 
    plot( ts, Rs[2],'r:' );  plot( ts, Os[2], 'b:' );   
    plot( ts, Fs[1],'r-' );  plot( ts, Fs[0], 'b.-' );   
    plot( ts, Fs[2],'r--' );   # G
    plot( ts, Fs[3],'r.-' );   # FTR
    plot( ts, Fs[4],'k.-' );   # FT 
    grid()

def map2D( X, U1, U2, f1, f2, n1, n2 ):
    #print " X: ",X
    M = zeros((2*n1+1,2*n1+1))
    for i in range(-n1,n1+1):
        d1 = array(U1)*(i*f1/n1)
        for j in range(-n2,n2+1):
            d2 = array(U2)*(j*f2/n2)
            M[i+n1,j+n2] = evalFitness( array(X)+d1 +d2 )
    return M

def plotMaps(irow,nrow, Gen):
    for i in range(nnodes):
        U1 = zeros(2*nnodes); U1[i       ]=1.0
        U2 = zeros(2*nnodes); U2[i+nnodes]=1.0
        #print " maping node",i," U1: ",U1," U2: ", U2
        subplot(nrow, nnodes, nnodes*irow+i+1 )
        mapa = map2D( Gen, U1, U2, 0.1, 0.1, 5, 5 )
        imshow(mapa, interpolation='bicubic', cmap='jet'); colorbar( )
        CS = contour(mapa, colors="g"); clabel(CS, inline=0.5, fontsize=8)

    
# ================ MAIN PROGRAM BODY =========================
figure(num=None, figsize=(20, 10))

GenHistory = []

#Gen = array([ A0es[0], A1es[0], A0es[1], A1es[1] ])
#Gen = MCrun( Gen, 10.0,  0.001, 100, 10  )

print " Initial Gen : ", Gen
Gen0 = array(Gen).copy()

Os,Rs,Fs = evalGen(Gen)

print "shape(ts)", shape(ts)
print "shape(ts)", shape(Fs[4])
#fitness = fitnesFunc(Fs)
plotTrj( Os,Rs,Fs, 0 )

#Simp = Simplex(evalFitness, Gen, [1, 1, 1, 1])

steps = ones(nnodes*2)*0.05
Simp = Simplex(evalFitness, Gen, steps )
#values, err, niter = SimplexOpt.minimize()
old_low = 10000000000
lastImprovement = 0
for i in range(0, 10000):
    converged, err,low,hi = Simp.simplexStep( 0.00001 )
    if converged:
        print " converged in ",i," steps " 
        break;
    if(low < old_low):
        lastImprovement = i
        old_low = low 
        subplot(2,5,5); plot( ts, Fs[4], '-', lw=0.25 );  grid()
        GenHistory.append(list(Simp.simplex[Simp.lowest]))
        print " new_low :  ", low, " iter: ", i, " err ", err
    if(i-lastImprovement)>(nnodes*8):
        print " Not able to improve => Exiting .... " 
        break;
print Simp.simplex[Simp.lowest]
subplot(2,5,5); plot( ts, Fs[4], 'r.-', lw=2 );  grid()

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

print " Initial gen ", Gen0
print " final   gen ", Gen

figure(num=None, figsize=(18, 4))
plotMaps(0,2, Gen0)
plotMaps(1,2, Gen )


show()