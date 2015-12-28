
from pylab import *
import matplotlib.ticker as ticker
from scipy.optimize import minimize,bisect
from PyKEP import lambert_problem
ax = subplot(111)
ax.xaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.xaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )
ax.yaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.yaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )

# ========= Problem definition

amax = 30
rmax = 2.0;   rmin=0.0001

phi = pi/2   
dt =0.5
R1 = array([1.0,0.0,0.0])

R20 = array([cos(phi),sin(phi),0.0])


# =========== Lambert function

def getLambert(r):
	R2 = R1 + R20*r
	a = NAN
	try:
		lambert_result =lambert_problem(R1,R2,dt )
		if lambert_result.is_reliable():
			a = lambert_result.get_a()[0]
	except:
		a= NAN
	#print " r,a ", r,a 
	return a
	

# ============ Pre optimize

x0 = bisect( getLambert, rmin,rmax, xtol=0.00001);  #   print " x0 = ", x0
print " bisect DONE "

def fitFractionalFunc( x0, A, B):
	a = (A[1]-B[1])/( 1/(x0-A[0]) - 1/(x0-B[0])  )
	b = A[1]- a/(x0-A[0])
	return a,b

ysc,y0 = fitFractionalFunc( x0, (rmin,getLambert(rmin)) , (x0*0.9,getLambert(x0*0.9)) );



'''
# NOT  GOOD
Gen = (y0, ysc, 0, 0, 0, 0, 0); print " Gen0 = ",Gen
def fitfuc(x, x0=0.55, P=( 0.3, 0.11, 0, 0, 0 ) ):
	ix1 = 1.0/(x0-x)
	ix2=ix1*ix1
	ix3=ix2*ix1
	ix4=ix3*ix1
	return P[0]+ P[1]*ix1 + P[2]*ix2 + P[3]*ix3 # + P[4]*ix4
'''

Gen = ( ysc + x0*y0, -y0, 0, 0, 0,0, 0, 0, 0, 0 ); print " Gen0 = ",Gen
def fitfuc(xx, x0=0.55, P=( 0.3, 0.11, 0, 0, 0, 0, 0, 0, 0, 0  ) ):
	x = (x0-xx)
	x2=x*x; x3=x2*x; 	x4=x3*x; x5=x4*x;
	p = P[0] + P[1]*x + P[2]*x2 + P[5]*x3 + P[7]*x4 + P[9]*x5
	q = x + P[3]*x2 + P[4]*x3 + P[6]*x4 + P[8]*x5
	return p/q


# =========== evaluate reference

print "=========== evaluate reference"

rs = arange(rmin,1.000001,0.05); nr=len(rs);
#rs = (rs**0.5)*x0*0.95
rs = rs*x0*0.95
 
ar = zeros(nr)
for ir in range(nr):
	ar[ir] = getLambert(rs[ir])

afit = fitfuc( rs, x0=x0, P=Gen)
plot(rs,afit)

plot(rs,ar)

# ========== Optimize ============

aref = ar.copy()
aref[aref>amax] = NaN
aref[aref<0]    = NaN



def evalFitness(Gen):
	afit = fitfuc( rs, x0=x0, P=Gen )
	err = nansum( ( afit - aref)**2 )
	#print err
	return err

res = minimize(evalFitness, Gen, method='BFGS', tol=0.0001, options={ 'maxiter':40, 'disp': False})
Gen = res.x
print " GenOpt = ", Gen

afit = fitfuc( rs, x0=x0, P=Gen)

plot(rs,aref,'o')
plot(rs,afit, lw=2)
# =========== Plotting

err = (afit-ar)
plot(rs,err*100)

ylim(-aref[-1]*1.2,aref[-1]*1.2)
axhline(0,ls='--')

figure()
plot(rs,abs(err)); 
yscale("log")

show()


