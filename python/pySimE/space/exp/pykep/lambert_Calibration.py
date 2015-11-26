
from pylab import *

from PyKEP import lambert_problem


xs = arange(-2,2,0.05); nx=len(xs); 
ys = arange(-2,2,0.05); ny=len(ys); 
R1 = array([0.0,1.0,0.0])
R2 = array([0.0,0.0,0.0])
dt =0.7

axy = zeros((nx,ny))
for ix in range(nx):
	for iy in range(ny):
		R2[0] = xs[ix]; R2[1] = ys[iy];
		try:
			lambert_result =lambert_problem(R1,R2,dt )
		except:
			continue
		#lambert_result = lambert_problem([1,0,0],[0,1,0],5 * pi / 2. )
		if lambert_result.is_reliable():
			a = lambert_result.get_a()[0]
			axy[ix,iy]= a
			print ix,iy," R2= ",R2," a= ", a 
		else:
			print ix,iy," R2= ",R2," failed "


imshow( axy, vmax=2.5,vmin=0, interpolation='nearest' )
colorbar()

show()


'''
lambert_result =lambert_problem( [0,1,0] ,[-0.5,0,0],1 )
a = lambert_result.get_a()[0]
print "a = ", a
'''


'''
dt =0.7
R1 = array([1.0,0.0,0.0]);   R2 = array([0.0,0.0,0.0])
phis = arange(0,1,0.25)*2*pi; nphi=len(phis); 
rs = arange(0,1,0.002); nr=len(rs); 
arphi = zeros((nphi,nr))
arphi[:,:] = NAN

for iphi in range(nphi):
	ca = cos(phis[iphi])
	sa = sin(phis[iphi])
	for ir in range(nr):
		R2[0] = ca*rs[ir]; R2[1] = sa*rs[ir];
		#print R2
		try:
			lambert_result =lambert_problem(R1,R2,dt )
			if lambert_result.is_reliable():
				a = lambert_result.get_a()[0]
				arphi[iphi,ir]= a
		except:
			#print "problem in lambert"
			continue
		#lambert_result = lambert_problem([1,0,0],[0,1,0],5 * pi / 2. )
	plot(rs,arphi[iphi])
	print iphi, " phi= ", phis[iphi]

ylim(-10,10)
#yscale('log')
axhline(0,ls='--')
'''


show()

