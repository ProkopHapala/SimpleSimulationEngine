
from PyKEP import *
from pylab import *

maxlines=1000000

asteroid_density = 2.0 # [ t/m^3 ]
min_radius = 500     # [ m ]  
max_radius = 10000   # [ m ]  


iline = 0
infile  = open('/home/prokop/Desktop/MPCORB.DAT')
#infile  = open('D:\\know\\Technika\\SpaceTech\\Asteroids\\Datafiles\\MPCORB.DAT')


for i in range(41):
	line = infile.readline()
	#print "---",i,"---",line

def toScale( x, xmin, xmax, nsteps ):
	if (x>=xmax or x<xmin):
		return -1
	else:
		return nsteps*(x-xmin)/(xmax-xmin)



Jupiter = planet_ss('Jupiter')
Mars    = planet_ss('Mars')
Earth   = planet_ss('Earth')

nsteps = 50
max_W  = (2*pi)*RAD2DEG
max_i  = (pi/8)*RAD2DEG

hist_all  = zeros((nsteps,nsteps*4))
hist_belt = zeros((nsteps,nsteps*4))
hist_NEA  = zeros((nsteps,nsteps*4))

data=[]

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	pl = planet_mpcorb(line)
	data.append([pl.radius,pl.orbital_elements[3],pl.orbital_elements[2]])
	ix = toScale(pl.orbital_elements[3]*RAD2DEG,0,max_W,nsteps*4 )
	iy = toScale(pl.orbital_elements[2]*RAD2DEG,0,max_i,nsteps )
	#print ix,iy
	if ( (ix>0) and (iy>0) and (pl.radius<max_radius) and (pl.radius>min_radius) and (pl.orbital_elements[0]<Jupiter.orbital_elements[0])  ):
		#data.append([pl.radius,pl.orbital_elements[3],pl.orbital_elements[2]])
		#print pl.radius,pl.orbital_elements[3],pl.orbital_elements[2]
		mass = asteroid_density * 4.18879020479 * pl.radius**3
		#mass = 1
		hist_all[iy,ix]+=mass; 
infile.close()

'''
data=transpose(array(data))
#scatter( data[1], data[2], s=log10(data[0]) )
scatter( data[1], data[2], s=1.0, marker='.' ,  alpha=0.25  )
xlim(2,3.5)
'''
figure(figsize=(15,5))

imshow( hist_all, extent=(0,max_W,0,max_i), origin='image', interpolation='nearest', aspect='equal' );
colorbar(); xlabel('semimajor axis [AU]'); ylabel('inclination [rad]'); title('all');

show()
