#!/usr/bin/env python

from pylab import *
from Poly6th_numeric import *


def CircularPV( R0, t ):
    cdO0  = sqrt(1/R0**3)
    return array([  R0 ,  cdO0*t  ]), array( [ 0.0, cdO0 ] )

T1= 4.0
T2= 4.0
Rinit = 1.0
Rend  = 1.5
Rmid  = sqrt( Rinit * Rend )

P0,V0 = CircularPV( Rinit , 0,   )
P1,V1 = CircularPV( Rmid , T1,   )
P2,V2 = CircularPV( Rend , T1+T2 )


print " P0,V0 ",P0,V0
print " P1,V1 ",P1,V1
print " P2,V2 ",P2,V2

#A0es  ,A1es_a = ForceEstimate_Cos ( P0, V0*T1, P1, V1*T1 )
#A1es_b,A2es   = ForceEstimate_Cos ( P1, V1*T2, P2, V2*T2 )

A0es  ,A1es_a = ForceEstimate_Cos ( P0, V0*T1, P1, V1*T1 )
A1es_b,A2es   = ForceEstimate_Cos ( P1, V1*T2, P2, V2*T2 )

#A0es  ,A1es_a = ForceEstimate_Linear ( P0, V0*T1, P1, V1*T1 )
#A1es_b,A2es   = ForceEstimate_Linear ( P1, V1*T2, P2, V2*T2 )

#A1es = (A1es_a + A1es_b) / 2
A1es = A1es_b
print " A0es ",A0es
print " A1es ",A1es, A1es_a, A1es_b
print " A2es ",A2es

# Evaluate Arrays
ts1 = arange(0,1.0001,0.05)
ts2 = arange(0,1.0001,0.05)


# Gen structure:
#     0    1    2     3    4   5    6    7    8    
#  [ PO1, VO1, VR1, AO0, AR0, AO1, AR1, AO2, AR2 ]

def evalGen ( Gen ):
    #print len(Gen)
    CRs1 = poly6Coefs(  P0[0],  V0[0]*T1 , Gen[4],  P1[0], Gen[2]*T1, Gen[6] )
    COs1 = poly6Coefs(  P0[1],  V0[1]*T1 , Gen[3], Gen[0], Gen[1]*T1, Gen[5] ) 
    CRs2 = poly6Coefs(  P1[0], Gen[2]*T2 , Gen[6],  P2[0],  V2[0]*T2, Gen[8] )
    COs2 = poly6Coefs( Gen[0], Gen[1]*T2 , Gen[5],  P2[1],  V2[1]*T2, Gen[7] ) 
    #CRs1 = poly6Coefs(  P0[0],  V0[0]*T1 , Gen[3]*T1**2,  P1[0], Gen[2]*T1, Gen[6]*T1**2 )
    #COs1 = poly6Coefs(  P0[1],  V0[1]*T1 , Gen[4]*T1**2, Gen[0], Gen[1]*T1, Gen[5]*T1**2 ) 
    #CRs2 = poly6Coefs(  P1[0], Gen[2]*T2 , Gen[6]*T2**2,  P2[0],  V2[0]*T2, Gen[8]*T2**2 )
    #COs2 = poly6Coefs( Gen[0], Gen[1]*T2 , Gen[5]*T2**2,  P2[1],  V2[1]*T2, Gen[7]*T2**2 ) 
    Os1,Rs1,Fs1 = evalPolarForces(ts1, T1,  CRs1, COs1)
    Os2,Rs2,Fs2 = evalPolarForces(ts2, T2,  CRs2, COs2)
    ts = concatenate( (ts1*T1, ts2[1:]*T2 + T1) )
    #print shape(Os1)
    Os = concatenate( (array(Os1), array(Os2)[:,1:]) ,axis=1 )
    Rs = concatenate( (array(Rs1), array(Rs2)[:,1:]) ,axis=1 )
    Fs = concatenate( (array(Fs1), array(Fs2)[:,1:]) ,axis=1 )  
    return ts,Os,Rs,Fs

def fitnesFunc( Fs ):
     return -sqrt( (( Fs[4][1:] + Fs[4][:-1]) * (ts[1:] - ts[:-1])).sum() )
    
stepSize = 1.0

def TryNew( GenBest, fitnessBest ):
    hit = False
    GenNew     =  GenBest[:] + (rand(len(GenBest))[:]-0.5)*stepSize 
    ts,Os,Rs,Fs   =  evalGen ( GenNew )
    fitnessNew = fitnesFunc(Fs)
    if(fitnessNew > fitnessBest ):
        hit = True
        GenBest     = GenNew
        fitnessBest = fitnessNew
        #print " Better is ",GenBest," fitness =  ",fitnessBest,
        #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
        subplot(2,5,5); plot( ts, Fs[4], '-', lw=0.25 );  grid()
    return GenBest, fitnessBest,hit

def plotTrj( ts, Os,Rs,Fs, i ):
    #subplot(2,5,1+5*i); plot( Os[0], Rs[0], '.-' );  grid()
    n2 = len(ts)/2
    subplot(2,5,1+5*i, polar=True); plot( Os[0][:n2], Rs[0][:n2], '.-g'); plot( Os[0][n2:], Rs[0][n2:], '.-k');  
    subplot(2,5,2+5*i);  plot( ts, Rs[0],'r-' ); plot( ts, Os[0], 'b-' );grid()
    subplot(2,5,3+5*i);  plot( ts, Rs[1],'r-' ); plot( ts, Os[1], 'b-' );grid() 
    subplot(2,5,4+5*i); 
    plot( ts, Rs[2],'r:' , lw=2);  plot( ts, Os[2], 'b:',lw=2 );   
    plot( ts, Fs[1],'r-' );  plot( ts, Fs[0], 'b.-' );   
    plot( ts, Fs[2],'r--' );   # G
    plot( ts, Fs[3],'r.-' );   # FTR
    plot( ts, Fs[4],'k.-' );   # FT 
    grid()


# ================ MAIN PROGRAM BODY =========================
figure(num=None, figsize=(18, 8))

fittness = -100000000

GenHistory = []

# Gen structure:
#     0    1    2     3    4   5    6    7    8    
#  [ PO1, VO1, VR1, AO0, AR0, AO1, AR1, AO2, AR2 ]

Gen = array([ P1[1], V1[1], V1[0], A0es[1], A0es[0], A1es[1], A1es[0], A2es[1], A2es[0] ])

ts,Os,Rs,Fs = evalGen(Gen)
fitness = fitnesFunc(Fs)
plotTrj( ts, Os,Rs,Fs, 0 )
print " Initail Gen ",Gen," fitness =  ",fitness
badInRow  = 0
fromMajor = 0
fitnessMajor = fitness
for i in range(1000):
    Gen,fitness,hit = TryNew( Gen, fitness )
    if(hit):
        #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
        print " fitness: ",fitness," step ",stepSize, " i: ", i, " bad ",badInRow," fromMajor",fromMajor
        GenHistory.append(Gen)
        badInRow = 0
        if(fitness-fitnessMajor)>0.025:
            fitnessMajor = fitness
            fromMajor    = 0
    if badInRow>30:
        stepSize *= 0.5
        badInRow = 0
    if fromMajor>100:
        break
    badInRow  += 1
    fromMajor += 1
ts, Os,Rs,Fs = evalGen(Gen)
plotTrj( ts, Os,Rs,Fs, 1 )

if len(GenHistory)>2:
    GenHistory = transpose(array(GenHistory ))
    subplot(2,5,10); 
    plot( GenHistory[0], '.-' ); 
    plot( GenHistory[1], '.-' ); 
    plot( GenHistory[2], '.-' ); 
    plot( GenHistory[3], '.-' ); 
else: 
    print "  NO BETTER SOLUTION FOUND "

show()