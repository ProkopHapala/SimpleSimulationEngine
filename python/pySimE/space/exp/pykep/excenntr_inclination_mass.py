
from PyKEP import *
from pylab import *

maxlines=1000000


asteroid_density = 2.0 # [ t/m^3 ]
max_radius = 10000     # [ m ]  

iline = 0
#infile  = open('/home/prokop/Desktop/MPCORB.DAT')
infile  = open('D:\\know\\Technika\\SpaceTech\\Asteroids\\Datafiles\\MPCORB.DAT')


for i in range(41):
	line = infile.readline()
	#print "---",i,"---",line

def toScale( x, xmin, xmax, nsteps ):
	if (x>=xmax or x<xmin):
		return -1
	else:
		return nsteps*(x-xmin)/(xmax-xmin)

nsteps = 10
max_e = pi/8
max_i = pi/8

hist_all = zeros((nsteps,nsteps))
hist_belt = zeros((nsteps,nsteps))
hist_NEA = zeros((nsteps,nsteps))

Jupiter = planet_ss('Jupiter')
Mars    = planet_ss('Mars')
Earth   = planet_ss('Earth')

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	pl = planet_mpcorb(line)
	ix = toScale(pl.orbital_elements[1],0,max_e,nsteps )
	iy = toScale(pl.orbital_elements[2],0,max_i,nsteps )
	if ( (ix>0) and(iy>0) and (pl.radius<max_radius)):
		mass = asteroid_density * 4.18879020479 * pl.radius**3
		hist_all[iy,ix]+=mass; 
		if (pl.orbital_elements[0] < Jupiter.orbital_elements[0]): # closer than Jupiter?
			hist_belt[iy,ix]+=mass
		if (pl.orbital_elements[0] < Mars.orbital_elements[0]): # closer than Mars?
			hist_NEA[iy,ix]+=mass 
infile.close()

figure(figsize=(15,5))



subplot(1,3,1); imshow( hist_all, extent=(0,max_e,0,max_i), origin='image', interpolation='nearest', aspect='equal' );
colorbar(); xlabel('eccentricity [rad]'); ylabel('inclination [rad]'); title('all');

subplot(1,3,2); imshow( hist_belt, extent=(0,max_e,0,max_i), origin='image', interpolation='nearest', aspect='equal' );
colorbar(); xlabel('eccentricity [rad]'); ylabel('inclination [rad]'); title('belt');

subplot(1,3,3); imshow( hist_NEA, extent=(0,max_e,0,max_i), origin='image', interpolation='nearest', aspect='equal' );
colorbar(); xlabel('eccentricity [rad]'); ylabel('inclination [rad]'); title('NEA');


show()
