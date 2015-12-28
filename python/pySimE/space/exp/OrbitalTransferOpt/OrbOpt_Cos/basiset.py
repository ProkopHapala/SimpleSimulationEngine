from pylab import *
#import numexpr as ne

# ============= Cosinus BaisSet ====================

def evalDiCos(i, ts, T):
	f1  = pi*(i); 	      f2 = pi*(i+2) 
	cs1 = cos(f1*ts); 	cs2 = cos(f2*ts)
	ss1 = sin(f1*ts); 	ss2 = sin(f2*ts)
	Y   =        cs1 -     cs2
	dY  = (   -f1*ss1 +     f2*ss2   )/( T    )
	ddY = (-f1*f1*cs1 +     f2*f2*cs2)/( T**2 )
	return array([Y,dY,ddY])

def evalDiCosBasisset( m, n, T ):
	ts = arange(0,1.000001,1.0/n)
	B  = []
	for i in range(m):
		B.append( evalDiCos(i, ts, T) )
	B = array(B)
	#print shape(B)
	B = transpose(B, (1,0,2) )
	#print shape(B)
	return B

# a - vector of amplitude coefficients
# B - matrix of basis fuctions for the coefficnets at time samples 0..T

'''
def evalSeriesInBasi(a,B):
	Y   = dot(a,B[0])
	dY  = dot(a,B[1])
	ddY = dot(a,B[2]) 
	return array([Y,dY,ddY])
'''

def evalSeriesInBasi(a,B):
	result = empty( (B.shape[0], B.shape[2]) )
	for i in xrange(B.shape[0]):
		dot(a, B[i], result[i])
	return result

# ============= Polygon4 BaisSet ====================

def poly4coefs_x0x1v0v1(x0, x1, v0, v1):
    A = zeros(4)
    A[0] =    x0
    A[1] =                v0
    A[2] = -3*x0 +3*x1 -2*v0 - v1  
    A[3] = +2*x0 -2*x1 +  v0 + v1
    return A
    
def evalPoly4( t, A):
    Y   = A[0] + A[1]*t +    A[2]*t**2 +    A[3]*t**3
    dY  =        A[1]   +  2*A[2]*t    +  3*A[3]*t**2
    ddY =                  2*A[2]      +  6*A[3]*t
    return array([Y,dY,ddY])

# =========== Evaluate Trajectory in polar coordinates ====================

def evalTrajectoryPolar( Rt0, Ot0, Bs, Rc, Oc ):
	Rt = Rt0 +  evalSeriesInBasi(Rc,Bs)
	Ot = Ot0 +  evalSeriesInBasi(Oc,Bs)
	Ft = evalPolarForces( Rt, Ot )
	#print " evalTrajectoryPolar shape Ot,Rt,Ft", shape(Ot), shape(Rt), shape(Ft)
	return Ot, Rt, Ft

def evalPolarForces( R, O ):	
	# numexpr doest seem to help it takes 3,644 vs. 1.910 with pure numpy
	# FIXME : check the sign of  F_R and G acording to wiki
	G    = -1.0 / (R[0]**2)                     # Gravitational force
	F_O  = R[0] * O[2] + 2 * R[1] * O[1]       # Angular Kinematic Force = Angular engine thrust 
	F_R  = R[2]        -     R[0] * O[1]**2    # Radial Kinematic force
	FTR  = F_R - G                             # Radial engine thrust
	FT2  = F_O**2 + FTR**2                     # Square of Total engine Trust Force ( corespons to propelant consuption for power limited variable specific impulse engine)	
	#print " shape F_O,F_R,G,FTR, FT2 ", shape(F_O),shape(F_R),shape(G),shape(FTR), shape(FT2)
	return [F_O,F_R,G,FTR, FT2]
	#return array([F_O,F_R,G,FTR, FT2])


'''
def evalPolarForces( R, O ):	
	# numexpr doest seem to help it takes 3,644 vs. 1.910 with pure numpy
	# FIXME : check the sign of  F_R and G acording to wiki
	R0 = R[0]; O0 = O[0];
	R1 = R[1]; O1 = O[1];
	R2 = R[2]; O2 = O[2];
	G    = ne.evaluate('-1.0 / (R0**2)'                  )   # Gravitational force
	F_O  = ne.evaluate('R0 * O2 + 2 * R1 * O1'   )   # Angular Kinematic Force = Angular engine thrust 
	F_R  = ne.evaluate('R2        -   R0 * O1**2')   # Radial Kinematic force
	FTR  = ne.evaluate('F_R - G'                         )   # Radial engine thrust
	FT2  = ne.evaluate('F_O**2 + FTR**2'                 )   # Square of Total engine Trust Force ( corespons to propelant consuption for power limited variable specific impulse engine)	
	#print " shape F_O,F_R,G,FTR, FT2 ", shape(F_O),shape(F_R),shape(G),shape(FTR), shape(FT2)
	return [F_O,F_R,G,FTR, FT2]
	#return array([F_O,F_R,G,FTR, FT2])
'''
