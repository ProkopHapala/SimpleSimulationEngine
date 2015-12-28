
from pylab import *

from PyKEP import epoch, DAY2SEC, planet_ss, AU, MU_SUN, lambert_problem
from PyKEP.orbit_plots import plot_planet, plot_lambert


plEarth = planet_ss('earth'); 
plMars  = planet_ss('mars');  

def opt_dt(tt1,tt2):
	t1 = epoch(tt1)
	t2 = epoch(tt2)
	#print t1
	dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC
	
	rE, vE = plEarth.eph(t1); vE=array(vE)
	rM, vM = plMars .eph(t2); vM=array(vM)	
	
	l = lambert_problem(rE,rM,dt,MU_SUN)
	
	vEl = array(l.get_v1()[0]); dvE = (vEl - vE)
	vMl = array(l.get_v2()[0]); dvM = (vMl - vM) 
	dvMTot = linalg.norm(dvM); dvETot= linalg.norm(dvE)
	dvTot = dvMTot+dvETot
	print " t1 " ,tt1," t2 ", tt2," dt ",(tt2-tt1)," dv ", dvTot
	return vE, vM, vEl, vMl

tt1min = 5000
tt1max = tt1min+365*10
#tt1max = 600
dttmin = 100
dttmax = 400

step = 5
dvMapEarth = []
dvMapMars  = []
for tt1 in range(tt1min,tt1max,step):
	dvRowE = []
	dvRowM = []
	for dtt in range(dttmin,dttmax,step):
		vE, vM, vEl, vMl = opt_dt(tt1,tt1+dtt)
		dvRowE.append(linalg.norm(vE-vEl))
		dvRowM.append(linalg.norm(vM-vMl))
	dvMapEarth.append(dvRowE) 
	dvMapMars .append(dvRowM)
dvMapEarth = array( dvMapEarth )
dvMapMars  = array( dvMapMars  )

#print shape(dvMap)
#print dvMap

#figure(figsize=(10,8))
#figure(figsize=(18,3))
figure(figsize=(20,6))


xtcks = range(tt1min,tt1max,365); #print xtcks
xtck_labes = [ str(epoch(tt))[:12] for tt in xtcks   ]; #print xtck_labes

subplot (3,1,1)
title('Total Earth to Mars delta-v [km/s]')
FF = transpose((dvMapEarth+dvMapMars)/1000)
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=10 , origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel('Departure date')
ylabel('Time of flight [days]')

subplot (3,1,2)
title('Earth departure delta-v [km/s]')
FF = transpose(dvMapEarth/1000)
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=7 , origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel('Departure date')
ylabel('Time of flight [days]')

subplot (3,1,3)
title('Mars arrival delta-v [km/s]')
FF = transpose(dvMapMars/1000)
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=5 , origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel('Departure date')
ylabel('Time of flight [days]')

#cfig = contour( FF , extent=( tt1min,tt1max, dttmin, dttmax ), colors='black', levels=arange(0,10,0.5))
#clabel(cfig, inline=1, fontsize=10)



tight_layout()

savefig('porkchop_Earth_Mars.png', transparent=True, bbox_inches='tight', pad_inches=0)
plt.show()

