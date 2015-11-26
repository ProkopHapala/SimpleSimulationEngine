
from pylab import *
import matplotlib.ticker as ticker

from PyKEP import lambert_problem


ax = subplot(111)
ax.xaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.xaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )
ax.yaxis.set_major_locator( ticker.MaxNLocator(nbins=10) )
ax.yaxis.set_minor_locator( ticker.AutoMinorLocator(n=10) )


dt =0.7
R1 = array([1.0,0.0,0.0]);   
 
rs = arange(0,1,0.02); nr=len(rs); 
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


plot(rs,ar)

def fitfuc(x, x0=0.55, y0=0.11, yscale=0.3 ):
	dx = x0-x
	#return yscale/(dx - 0.3*dx**2) + y0
	#return yscale/dx**0.9 + y0
	#return sign(dx)*yscale/abs(dx)**0.92 + y0
	return sign(dx)*yscale*(1-x**2)/abs(dx) + y0

afit = fitfuc(rs)
plot(rs,afit)

plot(rs,100*(afit-ar))

ylim(-10,10)
#yscale('log')
axhline(0,ls='--')



show()

