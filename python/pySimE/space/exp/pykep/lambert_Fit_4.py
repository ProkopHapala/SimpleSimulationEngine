
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
dt =0.1
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

x0 = bisect( getLambert, rmin,rmax, xtol=0.0001);  #   print " x0 = ", x0
print " bisect DONE "

def fitFractionalFunc( x0, A, B):
	a = (A[1]-B[1])/( 1/(x0-A[0]) - 1/(x0-B[0])  )
	b = A[1]- a/(x0-A[0])
	return a,b

ysc,y0 = fitFractionalFunc( x0, (rmin,getLambert(rmin)) , (x0*0.9,getLambert(x0*0.9)) );

Gen = (x0, ysc + x0*y0, -y0, 0, 0, 0,0, 0, 0 ); print " Gen0 = ",Gen

def fitfuc(x, P=(0.55, 0.3, 0.11, 0, 0, 0, 0, 0, 0  ) ):
	x2=x*x; x3=x2*x; 	x4=x3*x; 
	p = P[1] + P[2]*x + P[3]*x2 + P[6]*x3 + P[8]*x4
	q = P[0] -      x + P[4]*x2 + P[5]*x3 + P[7]*x4
	return p/q

# =========== evaluate reference

print "=========== evaluate reference"

rs = arange(rmin,1.000001,0.05); nr=len(rs);
#rs = (rs**0.5)*x0*0.95
rs = rs*x0*0.95
 
ar = zeros(nr)
for ir in range(nr):
	ar[ir] = getLambert(rs[ir])

afit = fitfuc( rs, Gen)
plot(rs,afit)

plot(rs,ar)

# ========== Optimize ============

aref = ar.copy()
aref[aref>amax] = NaN
aref[aref<0]    = NaN



def evalFitness(Gen):
	afit = fitfuc( rs, P=Gen )
	err = nansum( ( afit - aref)**2 )
	#print err
	return err

res = minimize(evalFitness, Gen, method='BFGS', tol=0.0001, options={ 'maxiter':60, 'disp': False})
Gen = res.x
print " GenOpt = ", Gen

afit = fitfuc(rs, Gen)

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


