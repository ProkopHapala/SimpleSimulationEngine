



from pylab import *

from PyKEP import epoch, DAY2SEC, planet_ss, AU, MU_SUN, lambert_problem
from PyKEP.orbit_plots import plot_planet, plot_lambert

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(0,0,0, color='y')

plEarth = planet_ss('earth'); 
plMars  = planet_ss('mars');  

def opt_dt(tt1,tt2):
	t1 = epoch(tt1)
	t2 = epoch(tt2)
	dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC
	
	rE, vE = plEarth.eph(t1); vE=array(vE)
	rM, vM = plMars .eph(t2); vM=array(vM)	
	
	l = lambert_problem(rE,rM,dt,MU_SUN)
	
	vEl = array(l.get_v1()[0]); dvE = (vEl - vE)
	vMl = array(l.get_v2()[0]); dvM = (vMl - vM) 
	'''
	print ""
	print " detlal-V at Earth: "
	print " Earth: ",vE 
	print " Ship:  ",vEl
	print " delta  ", dvE, linalg.norm(dvE)
	print ""
	print " detlal-V at Mars : "
	print " Mars:  ",  vM 
	print " Ship:  ",  vMl
	dvM = (vM - vMl)
	print " delta  ", dvM, linalg.norm(dvM)
	print " total delta-v  ", linalg.norm(dvM)+linalg.norm(dvE)
	'''
	print " dt ",(tt2-tt1)," dv ", linalg.norm(dvM)+linalg.norm(dvE)
	plot_planet(ax,plMars, t0=t2, color=(0.8,0.8,1),  units = AU)
	plot_lambert(ax,l      , color=(1,0,0),  units = AU)

tt1 = 450

plot_planet(ax,plEarth, t0=epoch(tt1), color=(0.8,0.8,1), units = AU)
plot_planet(ax,plMars,  t0=epoch(tt1), color=(0.8,0.8,1),  units = AU)


for dtt in range(200,350,10):
	opt_dt(tt1,tt1+dtt)

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
								
axisEqual3D(ax)

plt.show()

