#!/usr/bin/env python

from pylab                import *
from basiset              import *
from Simplex_optimization import Simplex
from Random_optimization  import MCBias_Run,MCBias2_Run

nnodes   = 8
nsamp    = 64
    
Gen = [0.0]*2*nnodes

def R2omega(R):
    return sqrt(1.0/R**3)

T  = 5.0;   R0 = 1.0;   R1 = 0.2;  
v0=R2omega(sqrt(1.0/R0**3));   v1=sqrt(1.0/R1**3);
ph = 0.4*T*sqrt(v0**2+v1**2)
P0=array( [ R0 , 0  ] );        V0=array( [ 1.0, v0 ] ); 
P1=array( [ R1 , ph ] );        V1=array( [ 0, v1 ] ); 

Bs = evalDiCosBasisset( nnodes, nsamp,  T)


#scalesR = 1.0/array(range(1,nnodes+1))
scalesO = 10.0/array(range(1,nnodes+1))
scalesR = ones(nnodes)
#scalesO = ones(nnodes)

#scalesR = 1.0/array(range(1,nnodes+1))**2
#scalesO = 2.0/array(range(1,nnodes+1))**1.5
#scalesR = 1.0/array(range(1,nnodes+1))
#scalesO = 2.0/array(range(1,nnodes+1))

ts         = arange(0,1.000001,1.0/nsamp ) 
timescales = b = matrix([1.0, 1.0/T, 1.0/T**2] ).transpose();
Rt0 = array(multiply( timescales, evalPoly4( ts, poly4coefs_x0x1v0v1(P0[0], P1[0], V0[0]*T, V1[0]*T) ) ) )
Ot0 = array(multiply( timescales, evalPoly4( ts, poly4coefs_x0x1v0v1(P0[1], P1[1], V0[1]*T, V1[1]*T) ) ) )
ts  *= T

nEvaluations = 0

maxThrust = 2.0
def fitnesFunc( Fs ):
	global nEvaluations;
	nEvaluations +=1
	fsum = 0
	tsum = 0
	#print "len(Fs) ", len(Fs[4])," len(ts) ", len(ts)
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

def map1D( X, U, f, n ):
	M = zeros(2*n+1)
	for i in range(-n,n+1):
		d = array(U)*(i*f/n)
		M[i+n] = evalFitness( array(X)+d )
	return M

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

#Gen0 = [ -2.73325403e-01 ,  1.65400287e-01 ,  1.33554628e-01,  -8.44902438e-02,  -3.02544694e-02 ,  8.73492900e-02 , -4.78252654e-02  , 1.40012492e-02,  -1.87602635e-01,  -5.13258160e-03,   3.23920204e-02 , -4.18545079e-02,   7.94684839e-03 , -1.45004222e-04  , 4.00265813e-03 , -2.23496019e-03 ]
			
#Gen0 = [-0.56200879 , 0.85362106, -0.45522688,  0.10272901,  0.09894236, -0.10673499,  0.07026536, -0.01868364, -0.10466653,  0.10724102, -0.05727427,  0.0184174, -0.00144521,  0.01436489, -0.00957926,  0.01494106]	

#Gen0 =  [-0.48470625,  0.67382507, -0.2976267,   0.03945806,  0.07932987, -0.0528763,  0.02676187, -0.00296634, -0.12791386,  0.0787777,  -0.03275983,  0.00563085, -0.00123082,  0.01072641, -0.00701855,  0.00919509]		

#Gen0 = array(Gen0)
#Gen = Gen0.copy()
print " Initial Gen : ", Gen
evalFitness( Gen )
Gen0 = array(Gen).copy()



Opt = True

if Opt:
	nEvaluations=0
	#Gen  = MCBias2_Run( evalFitness, Gen, 0.5, 0.01, 4*4*nnodes, 2*nnodes, 1000,  GenHistory,  wStep=0.5, fBias = 3.0,  kBias0=0.8, kBiasf = 0.97    ) # good
	Gen  = MCBias2_Run( evalFitness, Gen, 0.5, 0.01, 4*4*nnodes, 2*nnodes, 1000,  GenHistory,  wStep=0.7, fBias = 2.0, kBias0=0.8, kBiasf = 0.95  ) # good
	GenRnd = array(Gen).copy()
	print "===== nEvaluations : ", nEvaluations

#	nEvaluations=0
#	for i in range(5):
#		steps = ones(nnodes*2)*0.02*(5-i)
#		Gen      = Simplex_Run(Gen,steps, GenHistory)
#	GenSimp  = array(Gen).copy()
#	print "===== nEvaluations : ", nEvaluations
	
	print Gen

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

#	print " ===== Simplex    Fittness ", evalFitness( GenSimp )
#	plotTrj( Os,Rs,Fs, 3, "k" )
#	subplot(2,5,5); plot( ts, Fs[4], 'k-', lw=2 );  grid()

print " ===== Initial    Fittness ", evalFitness( Gen0 )
plotTrj( Os,Rs,Fs, 1, "r" )
subplot(2,5,5); autoscale(False); plot( ts, Fs[4], 'r-', lw=2 );  grid(), title("propelant \n consumption");


subplot(2,5,4); yticks(arange(0,100)/10.0); title(" phi coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i ]=1.0
	M = map1D(Gen,U,0.1,3)
	plot(log10(M), label=str(i));
legend() 

subplot(2,5,9); yticks(arange(0,100)/10.0); title(" R coef stiffnes ");
for i in range(nnodes):
	U = zeros(2*nnodes); U[i+nnodes]=1.0
	M = map1D(Gen,U,0.1,3)
	plot(log10(M), label=str(i)); 
legend()

#savefig("plost.png", bbox_inches='tight')


'''
figure(num=None, figsize=(20, 5))
plotMaps(0,2, Gen0)
plotMaps(1,2, Gen )
savefig("valley.png", bbox_inches='tight')
'''

show()
