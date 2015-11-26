
from PyKEP import *
from pylab import *

maxlines=1000000

asteroid_density = 2.0 # [ t/m^3 ]
min_radius = 500     # [ m ]  
max_radius = 10000   # [ m ]  

ia = 5
ix = 5
iy = 5
fname='D:\\know\\Technika\\SpaceTech\\Asteroids\\Datafiles\\sorted\\tile_a%02i_Lx%02i_Ly%02i.dat' %(ia,ix,iy); print "file: ",fname
infile = open(fname, 'r' )
infile.readline()

data    = []
data2   = []
data3   = []
data4   = []
lines   = []
chosen = []

x0=0.1;
y0=-0.1;

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	pl = planet_mpcorb(line)
	#ci = cos([pl.radius,pl.orbital_elements[2]) 
	'''
	si = sin(pl.orbital_elements[2])
	cW = cos(pl.orbital_elements[3])
	sW = sin(pl.orbital_elements[3])
	x  = si*cW
	y  = si*sW
	'''
	a  = pl.orbital_elements[0]/AU
	si = sin(pl.orbital_elements[2])
	cW = cos(pl.orbital_elements[3])
	sW = sin(pl.orbital_elements[3])
	X  = si*cW
	Y  = si*sW
	cw = cos(pl.orbital_elements[4]+pl.orbital_elements[3])
	sw = sin(pl.orbital_elements[4]+pl.orbital_elements[3])
	x  = pl.orbital_elements[1]*cw
	y  = pl.orbital_elements[1]*sw 
	if ( (pl.radius<max_radius) and (pl.radius>min_radius)  ):
		lines.append(line)
		#if pl.orbital_elements[1]<0.02:
		if (sqrt((x-x0)**2 + (y-y0)**2)<0.02):
			#print len(lines)-1
			ii=len(lines)-1
			chosen.append(ii)
			data2.append([a,x,y,pl.radius])
			data4.append([a,X,Y,pl.radius])
			#print ii,pl.orbital_elements[1],pl.orbital_elements[4]
			#print ii,pl.orbital_elements
		data.append([a,x,y,pl.radius])
		data3.append([a,X,Y,pl.radius])
infile.close()

N = len(lines)
print " number of asteroids: ",N

print chosen

figure(figsize=(15,5))
data=transpose(array(data))
data2=transpose(array(data2))
data3=transpose(array(data3))
data4=transpose(array(data4))

subplot(1,3,1); scatter(  data[1], data[2], s=1.0, marker='.' ,  alpha=0.5  ); axis('equal')
scatter( data2[1], data2[2], s=10.0, marker='o'  );
subplot(1,3,2); scatter(  data3[1], data3[2], s=1.0, marker='.' ,  alpha=0.5  ); axis('equal')
scatter( data4[1], data4[2], s=10.0, marker='o'  );
#subplot(1,3,1); 
#plot( data [1], data [2], '.' ); 
#plot( data2[1], data2[2], 'o' ); axis('equal')
#subplot(1,3,2); scatter( data[0], data[1], s=1.0, marker='.' ,  alpha=0.5  )
#subplot(1,3,3); scatter( data[0], data[2], s=1.0, marker='.' ,  alpha=0.5  )


N = len(lines)
print " number of asteroids: ",N

#ilist = [2,7,10,13]
#ilist = [12,23,29,41]
#ilist = [158,196,308]
ilist  = chosen

from PyKEP.orbit_plots import plot_planet, plot_lambert
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')


t1 = epoch(0)
n=len(ilist)
for i in range(n):
	pl = planet_mpcorb(lines[ilist[i]])
	#print ilist[i], pl.orbital_elements[1], pl.orbital_elements[4]
	#print ilist[i], pl.orbital_elements
	print "=====", lines[ilist[i]]
	plot_planet(ax,pl, t0=t1, color=(1.0-i/float(n),0.5,i/float(n)), legend=True, units = AU)
	#plot_planet(ax,pl, t0=t1, color=(1.0-i/float(n),0.5,i/float(n)), legend=False, units = AU)

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
								
axisEqual3D(ax)

show()
