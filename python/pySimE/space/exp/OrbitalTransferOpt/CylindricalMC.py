#!/usr/bin/env python

from pylab import *
from Poly6th_numeric import *

from Simplex2 import Simplex



T= 2.0

def CircularP0V0( R0, T, phi0, vr ):
    cdO0  = sqrt(1/R0**3)
    return array([  R0 ,  phi0  ]), array( [ vr, cdO0 ] )*T
    
P0,V0 = CircularP0V0( 0.2, T, 0.0  , 0 )
P1,V1 = CircularP0V0( 1.0, T, 3*pi , -1 )

'''
def CircularP0V0( R0,T, t ):
    cdO0  = sqrt(1/R0**3)
    return array([  R0 ,  cdO0*t  ]), array( [ 0.0, cdO0 ] )*T
P0,V0 = CircularP0V0( 0.5, T, 0.0 )
P1,V1 = CircularP0V0( 1.0, T, T*5 )
'''



'''
T= 3.0
P0,V0 = CircularP0V0( 1.0, T, 0.0 )
P1,V1 = CircularP0V0( 1.0, T, T )
'''
# A0es,A1es = ForceEstimate_Cos       ( P0, V0, P1, V1 )
A0es,A1es = ForceEstimate_Linear       ( P0, V0, P1, V1 )
print " A0es :", A0es, " A1es :", A1es

A0 = A0es
A1 = A1es

# Evaluate Arrays
ts = arange(0,1.0001,0.025)

def evalGen ( Gen ):
    CRs = poly6Coefs( P0[0], V0[0],  Gen[0],   P1[0], V1[0], Gen[1] )
    COs = poly6Coefs( P0[1], V0[1],  Gen[2],   P1[1], V1[1], Gen[3] ) 
    return evalPolarForces(ts, T,  CRs, COs)


maxThrust = 2.0
def fitnesFunc( Fs ):
    fsum  = 0
    tsum = 0
    for i in range(len(Fs[4])-1):
        dt=(ts[i+1]-ts[i])
        df=0.5*(Fs[4][i+1]+Fs[4][i])
        fsum+=df*dt
        tsum+=dt
        df_over = df-maxThrust
        #if(df_over>0):
        #    fsum+= (df_over**2) * dt   # penalty for overloading engine
    return -fsum
    #return -T*  sqrt((Fs[4]**2).sum()) /len(ts)
    

def TryNew( GenBest, fitnessBest, stepSize ):
    hit = False
    GenNew     =  GenBest[:] + (rand(4)[:]-0.5)*stepSize 
    Os,Rs,Fs   =  evalGen ( GenNew )
    fitnessNew = fitnesFunc(Fs)
    if(fitnessNew > fitnessBest ):
        hit = True
        GenBest     = GenNew
        fitnessBest = fitnessNew
        #print " Better is ",GenBest," fitness =  ",fitnessBest,
        #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
        subplot(2,5,5); plot( ts, Fs[4], '-', lw=0.25 );  grid()
    return GenBest, fitnessBest,hit

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

def MCrun( Gen, stepSize, accuracy, stopAfter, missRatio ):
    fitness = -100000000
    Os,Rs,Fs = evalGen(Gen)
    fitness = fitnesFunc(Fs)
    plotTrj( Os,Rs,Fs, 0 )
    print " Initail Gen ",Gen,
    print " fitness =  ",fitness
    badInRow  = 0
    fromMajor = 0
    fitnessMajor = fitness
    for i in range(10000):
        Gen,fitness,hit = TryNew( Gen, fitness, stepSize )
        if(hit):
            #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
            badInRow = 0
            stepSize *= 1.3
            if(fitness-fitnessMajor)>accuracy:
                print " fitness: ",fitness," step ",stepSize, " i: ", i, " bad ",badInRow," fromMajor",fromMajor
                GenHistory.append(Gen)
                fitnessMajor = fitness
                fromMajor    = 0
        if badInRow>missRatio:
            stepSize *= 0.5
            badInRow = 0
        if fromMajor>stopAfter:
            break
        badInRow  += 1
        fromMajor += 1
    return Gen



# ================ MAIN PROGRAM BODY =========================
figure(num=None, figsize=(20, 10))

GenHistory = []

Gen = array([ A0es[0], A1es[0], A0es[1], A1es[1] ])
Gen = MCrun( Gen, 10.0,  0.001, 100, 10  )

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