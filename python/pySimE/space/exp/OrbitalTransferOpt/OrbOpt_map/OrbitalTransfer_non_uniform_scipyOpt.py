#!/usr/bin/env python

from pylab import *

from Simplex_optimization import Simplex
from ClampedCubicSpline import *
from Random_optimization import MCBias_Run,MCBias2_Run
from scipy.optimize import minimize
import cProfile

nnodes   = 16
nsplines = nnodes + 1 
perNode  = 8
nsamp    = (nsplines)*perNode 
    
Gen= [0.0]*2*nnodes

def R2omega(R):
    return sqrt(1.0/R**3)

nEvaluations = 0

T  = 6.0
R0 = 1.0
R1 = 0.2
v0=R2omega(R0)
v1=R2omega(R1)

P0=array( [R0,0] );    V0=array( [-0.5,v0] ); 

vo = v0
ph = 0.0 

ti = zeros(nnodes+2)
for i in range(nsplines):
    tt=(i+1)/float(nsplines)
    ti[i+1]= 0.5*( 3*tt - tt**2 ) * T

print " len(ti) ",len(ti)

print " ti: ",ti

for i in range(nnodes):
    w=(i+1)/float(nsplines); mw= 1.0-w
    ri       = mw*R0 + w*R1
    Gen[i]   = ri
    vi       = R2omega(ri)
    dt       = ti[i+1]-ti[i]
    ph      += 0.5*(vo+vi)*dt
    Gen[i+nnodes] = ph
    vo       = vi

vi  = R2omega(R1)
dt  = ti[nnodes+1]-ti[nnodes]
ph += 0.5*(vo+vi)*dt
P1=array( [R1,ph] );    V1=array( [0.5,v1] ); 

print " P0 ", P0
print " P1 ", P1
print " V0 ", V0
print " V1 ", V1


ts = array( range(nsamp) )*T/float(nsamp)

def gen2controPoints( Gen ):
    n = (len(Gen)/2)
    ri = [0]*(n+2); oi = [0]*(n+2);
    #print " n ri oi ",n, len(oi), len(ri)
    ri[0]  =P0[0]; oi[0]  =P0[1];
    ri[n+1]=P1[0]; oi[n+1]=P1[1];
    for i in range(0,n):
        ri[i+1]=Gen[i  ]
        oi[i+1]=Gen[i+n]
    return ri,oi


def evalGen ( ti, Gen ):
    #print len(ti),len(ri)
    ri,oi =  gen2controPoints( Gen )
    #print " len ti,ri,oi: ",len(ti),len(ri),len(oi)
    SR, arn = Spline4_clamped( ti, ri, V0[0], V1[0])
    SO, aon = Spline4_clamped( ti, oi, V0[1], V1[1])
    return evalPolarForces( T, SR, SO, perNode)

maxThrust = 2.0
def fitnesFunc( Fs ):
    global nEvaluations;
    nEvaluations +=1
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
    ts,Os,Rs,Fs   =  evalGen ( ti, Gen )
    fitness =  fitnesFunc(Fs)
    return -fitness

def plotTrj( Os,Rs,Fs, i, clr="k" ):
    #subplot(2,5,1+5*i); plot( Os[0], Rs[0], '.-' );  grid()
    ri,oi= gen2controPoints( Gen)
    subplot(2,5,1, polar=True); plot( Os[0], Rs[0], '.-'+clr);  plot( oi, ri, 'o'+clr);  
    subplot(2,5,2); 
    print len(ti),len(ri),len(oi)
    plot( ts, Rs[0],'-'+clr ); plot( ts, Os[0], '--'+clr ); 
    plot( ti, ri,'o'+clr ); plot( ti, oi, 'x'+clr ); 
    grid()
    subplot(2,5,3);  plot( ts, Rs[1],'-'+clr ); plot( ts, Os[1], '--'+clr ); grid()
    subplot(2,5,5+i); 
    plot( ts, Rs[2],'r--' );  plot( ts, Os[2], 'b--' );   
    plot( ts, Fs[1],'r-' );   plot( ts, Fs[0],  'b-' );   
    plot( ts, Fs[2],'g-');   # G
    plot( ts, Fs[3],'m-');   # FTR
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
				
def map1D( X, U, f, n ):
	M = zeros(2*n+1)
	for i in range(-n,n+1):
		d = array(U)*(i*f/n)
		M[i+n] = evalFitness( array(X)+d )
	return M

def plotMaps(irow,nrow, Gen):
    for i in range(nnodes):
        U1 = zeros(2*nnodes); U1[i       ]=1.0
        U2 = zeros(2*nnodes); U2[i+nnodes]=1.0
        print " maping node",i," U1: ",U1," U2: ", U2
        subplot(nrow, nnodes, nnodes*irow+i+1 )
        mapa = map2D( Gen, U1, U2, 0.1, 0.1, 3, 3 )
        imshow(mapa, interpolation='bicubic', cmap='jet'); colorbar( )
        CS = contour(mapa, colors="g"); clabel(CS, inline=0.5, fontsize=8)


def TryNew( GenBest, fitnessBest, stepSize ):
    hit = False
    GenNew     =  GenBest[:] + (rand(nnodes*2)[:]-0.5)*stepSize 
    ts,Os,Rs,Fs   =  evalGen ( ti, GenNew )
    fitnessNew = fitnesFunc(Fs)
    #fitnessNew = evalFitness( GenNew )
    if(fitnessNew > fitnessBest ):
        hit = True
        GenBest     = GenNew
        fitnessBest = fitnessNew
        #print " Better is ",GenBest," fitness =  ",fitnessBest,
        #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
        subplot(2,5,5); plot( ts, Fs[4], '-', lw=0.25 );  grid()
    return GenBest, fitnessBest,hit

def Simplex_Run(Gen,steps, GenHistory):
    print 
    print " ========= Simplex Optimization ================= "
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
        if(i-lastImprovement)>(nnodes*16):
            print " Not able to improve => Exiting .... " 
            break;
    print Simp.simplex[Simp.lowest]
    return Simp.simplex[Simp.lowest]
    
# ================ MAIN PROGRAM BODY =========================
figure(num=None, figsize=(20, 10))

GenHistory = []


print " Initial Gen : ", Gen
ts, Os,Rs,Fs = evalGen( ti, Gen)
plotTrj( Os,Rs,Fs, 1, "r" )
Gen0 = array(Gen).copy()

#Gen    = MC_Run( Gen, 10.0,  0.01, 4*8*nnodes, 8*nnodes , GenHistory )
#Gen    = MC_Run( Gen, 10.0,  0.05, 8*4*nnodes, 4*nnodes , GenHistory )
#Gen  = MCBias_Run( evalFitness, Gen, 0.1, 0.001, 4*4*nnodes, 2*nnodes, 5000,  GenHistory )
#                                                                                            wStep=0.7, fBias = 3.0,  kBias0=0.5, kBiasf = 0.95
#Gen  = MCBias2_Run( evalFitness, Gen, 0.1, 0.0001, 4*4*nnodes, 2*nnodes, 5000,  GenHistory,     wStep=0.5, fBias = 1.5,  kBias0=0.8, kBiasf = 0.9     )
#Gen  = MCBias2_Run( evalFitness, Gen, 0.1, 0.001, 4*4*nnodes, 2*nnodes, 5000,  GenHistory,     wStep=0.7, fBias = 3.0,  kBias0=0.5, kBiasf = 0.95    ) # good


'''
nEvaluations=0
Gen  = MCBias2_Run( evalFitness, Gen, 0.4, 0.01, 4*4*nnodes, 2*nnodes, 5000,  GenHistory,     wStep=0.5, fBias = 3.0,  kBias0=0.6, kBiasf = 0.97    ) # good
print "===== nEvaluations : ", nEvaluations
GenRnd = array(Gen).copy()
ts, Os,Rs,Fs = evalGen( ti, Gen)
plotTrj( Os,Rs,Fs, 2, "g" )
nEvaluations=0

steps = ones(nnodes*2)*0.05
nEvaluations=0
Gen = Simplex_Run(Gen,steps, GenHistory)
print "===== nEvaluations : ", nEvaluations


ts, Os,Rs,Fs = evalGen( ti, Gen)
plotTrj( Os,Rs,Fs, 3, "k" )

if len(GenHistory)>2:
    GenHistory = transpose(array(GenHistory ))
    subplot(2,5,10);
    for i in range(nnodes):
        plot( GenHistory[i       ]-Gen0[i       ], 'r-' ); 
        plot( GenHistory[i+nnodes]-Gen0[i+nnodes], 'b-' );  
    #legend( bbox_to_anchor=(0.5, 1.00, 1., 0.000) )

ts, Os,Rs,Fs = evalGen( ti, Gen)
subplot(2,5,5); plot( ts, Fs[4], 'k-', lw=2 );  grid()

ts, Os,Rs,Fs = evalGen( ti, GenRnd)
subplot(2,5,5); plot( ts, Fs[4], 'g-', lw=2 );  grid()
'''

ts, Os,Rs,Fs = evalGen( ti, Gen0)
subplot(2,5,5); plot( ts, Fs[4], 'r-', lw=2 );  grid()


def minimizeRun():
	minimize(evalFitness, Gen, method='BFGS', tol=0.001, options={ 'maxiter':1000, 'disp': False})

nEvaluations=0
print " Minimizing ..... "
#res = minimize(evalFitness, Gen, method='Nelder-Mead', tol=0.01, options={ 'maxiter':1000, 'disp': True})
#res = minimize(evalFitness, Gen, method='Powell', tol=0.01, options={ 'maxiter':1000, 'disp': True})
res = minimize(evalFitness, Gen, method='BFGS', tol=0.001, options={ 'maxiter':1000, 'disp': False})
#cProfile.run('minimizeRun()',sort=1)
#res = minimize(evalFitness, Gen, method='COBYLA', tol=0.001, options={ 'maxiter':1000, 'disp': False})
#res = minimize(evalFitness, Gen, method='TNC', tol=0.001, options={ 'maxiter':1000, 'disp': False})
print res
Gen = res.x
print "===== nEvaluations : ", nEvaluations
ts, Os,Rs,Fs = evalGen( ti, Gen)
plotTrj( Os,Rs,Fs, 3, "k" )
subplot(2,5,5); plot( ts, Fs[4], 'g-', lw=2 );  grid()

print " Initial gen ", Gen0
print " final   gen ", Gen

savefig("plost.png", bbox_inches='tight')


'''
figure(num=None, figsize=(20, 5))
plotMaps(0,2, Gen0)
plotMaps(1,2, Gen )
savefig("valley.png", bbox_inches='tight')
'''


subplot(2,5,4); yticks(arange(0,100)/10.0); title(" R coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i ]=1.0
	M = map1D(Gen,U,0.2,10)
	plot(log10(M)); 

subplot(2,5,9); yticks(arange(0,1000)/100.0); title(" phi coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i+nnodes]=1.0
	M = map1D(Gen,U,0.2,10)
	plot(log10(M)); 

show()
