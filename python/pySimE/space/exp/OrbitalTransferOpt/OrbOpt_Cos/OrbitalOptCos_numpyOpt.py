#!/usr/bin/env python

from pylab                import *
from basiset              import *
from Simplex_optimization import Simplex
from Random_optimization  import MCBias_Run,MCBias2_Run
from scipy.optimize import minimize

import cProfile


nnodes   = 32
nsamp    = 256
    
Gen = [0.0]*2*nnodes

def R2omega(R):
    return sqrt(1.0/R**3)

#T  = 4.5;   R0 = 0.5;   R1 = 1.2;  
T  = 6.0;   R0 = 1.2;   R1 = 0.2;  
v0=R2omega(sqrt(1.0/R0**3));   v1=sqrt(1.0/R1**3);
ph = 0.4*T*sqrt(v0**2+v1**2)
P0=array( [ R0 , 0  ] );        V0=array( [  0.0, v0 ] ); 
P1=array( [ R1 , ph ] );        V1=array( [  0.0, v1 ] ); 

Bs = evalDiCosBasisset( nnodes, nsamp,  T)

scalesR = ones(nnodes)
scalesO = ones(nnodes)

#scalesR = 0.2*ones(nnodes)
#scalesO = 10.0/array(range(1,nnodes+1))

#scalesR = 1.0/array(range(1,nnodes+1))
#scalesO = 2.0/array(range(1,nnodes+1))

ts         = arange(0,1.000001,1.0/nsamp ) 
timescales = b = matrix([1.0, 1.0/T, 1.0/T**2] ).transpose();
Rt0 = array(multiply( timescales, evalPoly4( ts, poly4coefs_x0x1v0v1(P0[0], P1[0], V0[0]*T, V1[0]*T) ) ) )
Ot0 = array(multiply( timescales, evalPoly4( ts, poly4coefs_x0x1v0v1(P0[1], P1[1], V0[1]*T, V1[1]*T) ) ) )
ts  *= T
dts  = ts[1:] - ts[:-1] 

nEvaluations = 0

maxThrust = 2.0
def fitnesFunc( Fs ):
	global nEvaluations;
	nEvaluations +=1
	fsum = 0
	tsum = 0
	#print "len(Fs) ", len(Fs[4])," len(ts) ", len(ts)
	'''
	for i in range(len(Fs[4])-1):
		dt=(ts[i+1]-ts[i])
		df=0.5*(Fs[4][i+1]+Fs[4][i])
		fsum+=df*dt
		tsum+=dt
		#df_over = df-maxThrust
		#if(df_over>0):
		#    fsum+= (df_over**2) * dt   # penalty for overloading engine
	return -sqrt(fsum/tsum)
	#return -T*  sqrt((Fs[4]**2).sum()) /len(ts)
	'''
	#return -sqrt(((Fs[4,1:]+Fs[4,:-1])*dts).sum()/T) 
	return -sqrt(((Fs[4][1:]+Fs[4][:-1])*dts).sum()/T) 
    
def evalFitness( Gen ):
	global Os,Rs,Fs
	cR = Gen[nnodes:] * scalesR
	cO = Gen[:nnodes] * scalesO
	Os,Rs,Fs =  evalTrajectoryPolar( Rt0, Ot0, Bs, cR, cO )
	#print " evalFitness shape Os,Rs,Fs", shape(Rs),shape(Rs), shape(Fs)
	fitness  =  fitnesFunc(Fs)
	return -fitness

def plotTrj( Os,Rs,Fs, i, clr="k" ):
	print shape(ts), shape(Rs), shape(Rs), shape(Fs)
	subplot(2,5,1, polar=True); plot( Os[0], Rs[0], '-'+clr); title(" Trajectory ");
	subplot(2,5,2); 	plot( ts, Rs[0],'-'+clr ); plot( ts, Os[0], '--'+clr ); grid(); title(" Position ");
	subplot(2,5,3);  plot( ts, Rs[1],'-'+clr ); plot( ts, Os[1], '--'+clr ); grid(); title(" Velocity  ");
	subplot(2,5,5+i); 
	plot( ts, Rs[2],'r--' );  plot( ts, Os[2], 'b--' );   
	plot( ts, Fs[1],'r-' );   plot( ts, Fs[0], 'b-' );   
	plot( ts, Fs[2],'g-');   # G
	plot( ts, Fs[3],'m-');   # FTR
	plot( ts, sqrt(Fs[4]),'k.-' );   # FT 
	title(" acclereations ");
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
        print " maping node",i," U1: ",U1," U2: ", U2
        subplot(nrow, nnodes, nnodes*irow+i+1 )
        mapa = map2D( Gen, U1, U2, 0.1, 0.1, 3, 3 )
        imshow(mapa, interpolation='bicubic', cmap='jet'); colorbar( )
        CS = contour(mapa, colors="g"); clabel(CS, inline=0.5, fontsize=8)

def map1D( X, U, f, n ):
	M = zeros(2*n+1)
	for i in range(-n,n+1):
		d = array(U)*(i*f/n)
		M[i+n] = evalFitness( array(X)+d )
	return M

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
		converged, err,low,hi = Simp.simplexStep( 0.0001 )
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
evalFitness( Gen )
Gen0 = array(Gen).copy()

Opt = True

def minimizeRun():
	minimize(evalFitness, Gen, method='BFGS', tol=0.001, options={ 'maxiter':1000, 'disp': False})

if Opt:
	'''
	nEvaluations=0
	Gen  = MCBias2_Run( evalFitness, Gen, 0.5, 0.01, 4*4*nnodes, 2*nnodes, 2000,  GenHistory,  wStep=0.5, fBias = 3.0,  kBias0=0.8, kBiasf = 0.97    ) # good
	GenRnd = array(Gen).copy()
	print "===== nEvaluations : ", nEvaluations

	steps = ones(nnodes*2)*0.05
	nEvaluations=0
	Gen      = Simplex_Run(Gen,steps, GenHistory)
	GenSimp  = array(Gen).copy()
	print "===== nEvaluations : ", nEvaluations
	'''	


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
	evalFitness( Gen)
	plotTrj( Os,Rs,Fs, 3, "k" )
	subplot(2,5,5); plot( ts, Fs[4], 'k-', lw=2 );  grid()



'''
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
    10132   0.066    0.000    1.580    0.000 basiset.py:29(evalSeriesInBasi) # not optimized
    10132   0.103    0.000    0.514    0.000 basiset.py:36(evalSeriesInBasi) #     optimized

    5066    0.439    0.000    1.316    0.000 basiset.py:67(evalPolarForces)  # not optimized
    5032    0.350    0.000    0.350    0.000 basiset.py:67(evalPolarForces)  #     optimized ( no array() call )

    31201    0.425    0.000    0.425    0.000 {numpy.core._dotblas.dot}      #  not optimized
    31001    0.335    0.000    0.335    0.000 {numpy.core._dotblas.dot}      #      optimized
				
				
    10064    0.017    0.000    0.017    0.000 {numpy.core.multiarray.empty} #      optimized

    15200    1.971    0.000    1.971    0.000 {numpy.core.multiarray.array} #  not optimized				
    2        0.000    0.000    0.000    0.000 {numpy.core.multiarray.array} #      optimized
				
     5066    0.126    0.000    3.023    0.001 basiset.py:60(evalTrajectoryPolar)
     5066    0.111    0.000    1.941    0.000 basiset.py:60(evalTrajectoryPolar)
     5032    0.089    0.000    0.880    0.000 basiset.py:60(evalTrajectoryPolar)

     5066    0.083    0.000    3.278    0.001 OrbitalOptCos_numpyOpt.py:68(evalFitness) # evalPolarForces  optimized
     5066    0.083    0.000    2.197    0.000 OrbitalOptCos_numpyOpt.py:68(evalFitness) # evalSeriesInBasi optimized				
     5032    0.077    0.000    1.085    0.000 OrbitalOptCos_numpyOpt.py:69(evalFitness) # optimized both

  Without optimization
    15200    1.971    0.000    1.971    0.000 {numpy.core.multiarray.array}
     5066    0.445    0.000    1.317    0.000 basiset.py:67(evalPolarForces)
    31201    0.425    0.000    0.425    0.000 {numpy.core._dotblas.dot}
     5066    0.142    0.000    0.173    0.000 OrbitalOptCos_numpyOpt.py:48(fitnesFunc)
     5066    0.126    0.000    3.023    0.001 basiset.py:60(evalTrajectoryPolar)
      149    0.089    0.001    3.284    0.022 optimize.py:536(approx_fprime)
     5066    0.083    0.000    3.278    0.001 OrbitalOptCos_numpyOpt.py:68(evalFitness)
    10132    0.066    0.000    1.580    0.000 basiset.py:29(evalSeriesInBasi)
     5066    0.031    0.000    0.031    0.000 {method 'sum' of 'numpy.ndarray' objects}
  With optimization					
     5068    0.878    0.000    0.878    0.000 {numpy.core.multiarray.array}
     5066    0.439    0.000    1.316    0.000 basiset.py:67(evalPolarForces)
    31206    0.405    0.000    0.405    0.000 {numpy.core._dotblas.dot}
     5066    0.143    0.000    0.173    0.000 OrbitalOptCos_numpyOpt.py:48(fitnesFunc)
     5066    0.111    0.000    1.941    0.000 basiset.py:60(evalTrajectoryPolar)
    10132    0.103    0.000    0.514    0.000 basiset.py:36(evalSeriesInBasi)
      149    0.086    0.001    2.231    0.015 optimize.py:536(approx_fprime)
     5066    0.083    0.000    2.197    0.000 OrbitalOptCos_numpyOpt.py:68(evalFitness)
   More optimized
     5032    0.350    0.000    0.350    0.000 basiset.py:67(evalPolarForces)
    31001    0.335    0.000    0.335    0.000 {numpy.core._dotblas.dot}
     5032    0.104    0.000    0.128    0.000 OrbitalOptCos_numpyOpt.py:48(fitnesFunc)
    10064    0.098    0.000    0.441    0.000 basiset.py:36(evalSeriesInBasi)
     5032    0.089    0.000    0.880    0.000 basiset.py:60(evalTrajectoryPolar)
     5032    0.077    0.000    1.085    0.000 OrbitalOptCos_numpyOpt.py:69(evalFitness)
      148    0.074    0.001    1.138    0.008 optimize.py:536(approx_fprime)
        1    0.028    0.028    1.223    1.223 optimize.py:734(_minimize_bfgs)
     5032    0.024    0.000    0.024    0.000 {method 'sum' of 'numpy.ndarray' objects}
    10064    0.017    0.000    0.017    0.000 {numpy.core.multiarray.empty}
'''


'''
if len(GenHistory)>2:
	GenHistory = transpose(array(GenHistory ))
	subplot(2,5,10);
	for i in range(nnodes):
		plot( GenHistory[i       ]-Gen0[i       ], 'r-' ); 
		plot( GenHistory[i+nnodes]-Gen0[i+nnodes], 'b-' );  
		#legend( bbox_to_anchor=(0.5, 1.00, 1., 0.000) )

if Opt:
	print " ===== Random     Fittness ", evalFitness( GenRnd )
	plotTrj( Os,Rs,Fs, 2, "g" )
	subplot(2,5,5); plot( ts, Fs[4], 'g-', lw=2 );  grid()

	print " ===== Simplex    Fittness ", evalFitness( GenSimp )
	plotTrj( Os,Rs,Fs, 3, "k" )
	subplot(2,5,5); plot( ts, Fs[4], 'k-', lw=2 );  grid()
'''

print " ===== Initial    Fittness ", evalFitness( Gen0 )
plotTrj( Os,Rs,Fs, 1, "r" )
subplot(2,5,5); autoscale(False); plot( ts, Fs[4], 'r-', lw=2 );  grid(), title("propelant \n consumption");


'''
subplot(2,5,4); yticks(arange(0,100)/10.0); title(" phi coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i ]=1.0
	M = map1D(Gen,U,0.2,3)
	plot(log10(M)); 

subplot(2,5,9); yticks(arange(0,100)/10.0); title(" R coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i+nnodes]=1.0
	M = map1D(Gen,U,0.2,3)
	plot(log10(M)); 
	
savefig("plots_BFGS.png", bbox_inches='tight')
'''


'''
figure(num=None, figsize=(20, 5))
plotMaps(0,2, Gen0)
plotMaps(1,2, Gen )
savefig("valley.png", bbox_inches='tight')
'''

show()
