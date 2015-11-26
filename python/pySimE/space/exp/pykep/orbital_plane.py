
from PyKEP import *
from pylab import *

maxlines=1000000

asteroid_density = 2.0 # [ t/m^3 ]
min_radius = 500     # [ m ]  
max_radius = 10000   # [ m ]  


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



Jupiter = planet_ss('Jupiter')
Mars    = planet_ss('Mars')
Earth   = planet_ss('Earth')

nsteps = 10
max_X  = 0.3

hist_all  = zeros((nsteps,nsteps))
hist_belt = zeros((nsteps,nsteps))
hist_NEA  = zeros((nsteps,nsteps))

data=[]

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	pl = planet_mpcorb(line)
	#ci = cos([pl.radius,pl.orbital_elements[2]) 
	si = sin(pl.orbital_elements[2])
	cW = cos(pl.orbital_elements[3])
	sW = sin(pl.orbital_elements[3])
	x  = si*cW
	y  = si*sW
	#z  = ci
	#data.append([pl.radius,x,y])
	ix = toScale(x,-max_X,max_X,nsteps )
	iy = toScale(y,-max_X,max_X,nsteps )
	#print ix,iy
	if ( (ix>0) and (iy>0) and (pl.radius<max_radius) and (pl.radius>min_radius) and (pl.orbital_elements[0]<Jupiter.orbital_elements[0])  ):
		#data.append([pl.radius,pl.orbital_elements[3],pl.orbital_elements[2]])
		#print pl.radius,pl.orbital_elements[3],pl.orbital_elements[2]
		mass = asteroid_density * 4.18879020479 * pl.radius**3
		mass = 1
		hist_all[iy,ix]+=mass; 
infile.close()
'''
data=array(data)
print data
data=transpose(data)
#scatter( data[1], data[2], s=log10(data[0]) )
scatter( data[1], data[2], s=1.0, marker='.' ,  alpha=0.5  )
#scatter( data[1], data[2]  )
axis('equal')
'''
figure(figsize=(10,10))
imshow( hist_all, extent=(-max_X,max_X,-max_X,max_X), origin='image', interpolation='nearest', aspect='equal' );
colorbar(); xlabel('Lx'); ylabel('Ly]'); title('all');

show()
