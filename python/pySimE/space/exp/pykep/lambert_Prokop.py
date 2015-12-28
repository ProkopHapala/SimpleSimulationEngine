
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from PyKEP import epoch, DAY2SEC, planet_ss, AU, MU_SUN, lambert_problem
from PyKEP.orbit_plots import plot_planet, plot_lambert


mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')




t1 = epoch(0)
t2 = epoch(740)
dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

ax.scatter(0,0,0, color='y')

pl = planet_ss('earth')
plot_planet(ax,pl, t0=t1, color=(0.8,0.8,1), legend=True, units = AU)
rE,vE = pl.eph(t1)

pl = planet_ss('mars')
plot_planet(ax,pl, t0=t2, color=(0.8,0.8,1), legend=True, units = AU)
rM, vM = pl.eph(t2)

l = lambert_problem(rE,rM,dt,MU_SUN)

nmax = l.get_Nmax()
print "max number of revolutions",nmax

plot_lambert(ax,l      , color=(1,0,0), legend=True, units = AU)
for i in range(1,nmax*2+1):
	print i
	plot_lambert(ax,l,sol=i, color=(1,0,i/float(nmax*2)), legend=True, units = AU)



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

