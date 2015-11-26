
from PyKEP import *
from pylab import *

maxlines=1000000

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

data=zeros((3,maxlines))


nsteps = float(50)
max_log10_R  = float(5)

data_all  = zeros(nsteps)
data_NEA  = zeros(nsteps)
data_belt = zeros(nsteps)

Jupiter = planet_ss('Jupiter')
Mars    = planet_ss('Mars')
Earth   = planet_ss('Earth')

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	#print line
	pl = planet_mpcorb(line)
	ix = toScale( log10(pl.radius),0,max_log10_R,nsteps )
	#print pl.radius, log10(pl.radius),max_log10_R, ix
	if (ix>0):
		data_all[ix]+=1;
		if (pl.orbital_elements[0] < Jupiter.orbital_elements[0]): # closer than Jupiter?
			data_belt[ix]+=1
		if (pl.orbital_elements[0] < Mars.orbital_elements[0]): # closer than Mars?
			data_NEA[ix]+=1 
infile.close()

xs = arange(nsteps)*(max_log10_R/nsteps)
plot( xs, data_all  ,'.-', label='all')
plot( xs, data_belt ,'.-', label='belt')
plot( xs, data_NEA  ,'.-', label='NEA')

savetxt( 'size_histogram.dat', transpose(array([xs,data_all,data_belt,data_NEA ]) ) )



legend()
#ax.set_xscale('log')
yscale('log')

xlabel('log10 (radius) [m]')
ylabel('count')

show()
