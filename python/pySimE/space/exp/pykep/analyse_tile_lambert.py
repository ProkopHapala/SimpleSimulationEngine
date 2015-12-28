
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
infile  = open(fname, 'r' )
infile.readline()

data  = []
lines = [] 

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	pl = planet_mpcorb(line)
	#ci = cos([pl.radius,pl.orbital_elements[2]) 
	si = sin(pl.orbital_elements[2])
	cW = cos(pl.orbital_elements[3])
	sW = sin(pl.orbital_elements[3])
	x  = si*cW
	y  = si*sW
	a  = pl.orbital_elements[0]/AU
	if ( (pl.radius<max_radius) and (pl.radius>min_radius)  ):
		lines.append(line)
		data.append([a,x,y,pl.radius])
infile.close()

'''
figure(figsize=(15,5))
data=array(data)
print data
data=transpose(data)

subplot(1,3,1); scatter( data[1], data[2], s=1.0, marker='.' ,  alpha=0.5  ); axis('equal')
subplot(1,3,2); scatter( data[0], data[1], s=1.0, marker='.' ,  alpha=0.5  )
subplot(1,3,3); scatter( data[0], data[2], s=1.0, marker='.' ,  alpha=0.5  )
'''

N = len(lines)
print " number of asteroids: ",N


pl1= planet_mpcorb(lines[130] ); print pl1.name
pl2= planet_mpcorb(lines[170] ); print pl1.name

def opt_dt(tt1,tt2, pl1, pl2):
	t1 = epoch(tt1)
	t2 = epoch(tt2)
	dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC
	r1, v1 = pl1.eph(t1); v1=array(v1)
	r2, v2 = pl2.eph(t2); v2=array(v2)	
	l = lambert_problem(r1,r2,dt,MU_SUN)
	v1l = array(l.get_v1()[0]); dv1 = (v1l - v1)
	v2l = array(l.get_v2()[0]); dv2 = (v2l - v2) 
	dv1Tot = linalg.norm(dv1); dv2Tot= linalg.norm(dv2)
	dvTot = dv1Tot+dv2Tot
	print " t1 " ,tt1," t2 ", tt2," dt ",(tt2-tt1)," dv ", dvTot
	return v1, v2, v1l, v2l

#tt1min = 8000
tt1min = 365*100
#tt1min = 0
tt1max = tt1min+365*50
#tt1max = 600
dttmin = 100
dttmax = 3000

step = 50
dvMap1  = []
dvMap2 = []
for tt1 in range(tt1min,tt1max,step):
	dvRow1 = []
	dvRow2 = []
	for dtt in range(dttmin,dttmax,step):
		v1, v2, v1l, v2l = opt_dt(tt1,tt1+dtt, pl1, pl2)
		dvRow1.append(linalg.norm(v1-v1l))
		dvRow2.append(linalg.norm(v2-v2l))
	dvMap1.append(dvRow1) 
	dvMap2.append(dvRow2)
dvMap1 = array( dvMap1 )
dvMap2 = array( dvMap2  )


figure(figsize=(16,16))
xtcks = range(tt1min,tt1max,365); #print xtcks
#xtck_labes = [ str(epoch(tt))[:12] for tt in xtcks   ]; #print xtck_labes
xtck_labes = [ str(epoch(tt))[2:4] for tt in xtcks   ]; #print xtck_labes


dvmax=2

subplot (3,1,1)
title('Total delta-v [km/s]')
FF = transpose((dvMap1+dvMap2)/1000); vmin=FF.min()
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=vmin+dvmax , origin='image', interpolation='bicubic', aspect='equal');
#imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel(pl1.name+'Departure date')
ylabel('Time of flight [days]')

subplot (3,1,2)
title(pl1.name+' delta-v [km/s]')
FF = transpose(dvMap1/1000); vmin=FF.min()
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=vmin+dvmax , origin='image', interpolation='bicubic', aspect='equal');
#imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel('Departure date')
ylabel('Time of flight [days]')

subplot (3,1,3)
title(pl2.name+' delta-v [km/s]')
FF = transpose(dvMap2/1000); vmin=FF.min()
imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ), vmax=vmin+dvmax, origin='image', interpolation='bicubic', aspect='equal');
#imshow( FF, extent=( tt1min,tt1max, dttmin, dttmax ),  origin='image', interpolation='bicubic', aspect='equal');
colorbar(use_gridspec=True, shrink=0.9, pad = 0.005, fraction = 0.005 );
xticks( xtcks, xtck_labes  )
xlabel('Departure date')
ylabel('Time of flight [days]')

tight_layout()
#savefig('porkchop_Earth_Saturn.png', transparent=True, bbox_inches='tight', pad_inches=0)

show()
