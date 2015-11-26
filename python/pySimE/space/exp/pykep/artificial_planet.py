
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from PyKEP import *
from PyKEP.orbit_plots import plot_planet, plot_lambert

mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.scatter(0,0,0, color='y')

earth = planet(
epoch(54000,epoch.epoch_type.MJD),
(9.9998805e-01 * AU, 1.6716812e-02, 8.8543531e-04 * DEG2RAD, 1.7540648e+02 * DEG2RAD, 2.8761578e+02 * DEG2RAD, 2.5760684e+02 * DEG2RAD), 
MU_SUN, 398600.4418e9, 6378000, 6900000,  'Earth' )


pls = []
for i in range(5):
	myPlanet = planet(
	epoch(54000,epoch.epoch_type.MJD), (
	1.5 * AU,                    # a	semimajor axis
	0.5 ,               # e	eccentricity
	30 * DEG2RAD,     # i	inclination
	i*60 * DEG2RAD,     # W	longitude of accesing node
	0*60 * DEG2RAD,     # w	argument of periapsis
	0 * DEG2RAD      # M	mean anomaly
	), MU_SUN, 398600.4418e9, 6378000, 6900000,  'testPlanet' )
	pls.append(myPlanet)

t1 = epoch(0)
plot_planet(ax,earth,    t0=t1, color=(0.8,0.8,0.8), legend=True, units = AU)
for i in range(len(pls)):
	plot_planet(ax,pls[i], t0=t1, color=(1-i/5.0,0,i/5.0), legend=True, units = AU)



















# =================================================================================

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