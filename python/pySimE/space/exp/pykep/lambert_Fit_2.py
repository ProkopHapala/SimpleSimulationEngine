
from pylab import *
import matplotlib.ticker as ticker

from PyKEP import lambert_problem


ax = subplot(111)
ax.xaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.xaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )
ax.yaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.yaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )


amax = 10

dt =0.7
R1 = array([1.0,0.0,0.0]);   
 
rs = arange(0,1.0,0.02); nr=len(rs); 
rs = rs**0.5
ar = zeros(nr)
ar[:] = NAN

phi = 0.000000001
#phi = pi
phi = pi/2
R20 = array([cos(phi),sin(phi),0.0])

print R1
print R20

for ir in range(nr):
	R2 = R20*rs[ir];
	try:
		lambert_result =lambert_problem(R1,R2,dt )
		if lambert_result.is_reliable():
			a = lambert_result.get_a()[0]
			ar[ir]= a
	except:
		#print "problem in lambert"
		continue
		#lambert_result = lambert_problem([1,0,0],[0,1,0],5 * pi / 2. )


def fitfuc(x, P=(0.55, 0.11, 0.3, 0, 0, 0, 0  ) ):
	#return sign(dx)*yscale*(1-x**2)/abs(dx) + y0
	#return P[2]/( P[0] - x ) + P[1]
	return (P[2]+x*P[3]+x*x*P[5])/( P[0] - x + P[4]*x**2 + P[6]*x**3  ) + P[1]
	#return ( P[2]+ x*P[5] + P[3]*x**2 )/( P[0] - x + P[4]*x**2  ) + P[1]



aref = ar.copy()
aref[aref>amax] = NaN
aref[aref<0]    = NaN

Gen = (0.55, 0.11, 0.3, 0, 0, 0, 0 )
#Gen = (0.55, 0.11, 0.3, 0, 0, 0 )





# ========== Optimize ============
from scipy.optimize import minimize


def evalFitness(Gen):
	afit = fitfuc( rs, P=Gen )
	err = nansum( ( afit - aref)**2 )
	print err
	return err

res = minimize(evalFitness, Gen, method='BFGS', tol=0.001, options={ 'maxiter':100, 'disp': False})
Gen = res.x
print " Optimalized Gen ", Gen


afit = fitfuc(rs, Gen)

plot(rs,ar)
plot(rs,aref,'o')
plot(rs,afit, lw=2)
plot(rs,100*(afit-ar))

ylim(-10,30)
axhline(0,ls='--')


show()

