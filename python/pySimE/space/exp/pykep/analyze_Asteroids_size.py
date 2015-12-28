
from PyKEP import *
from pylab import *

maxlines=10000

iline = 0
infile  = open('/home/prokop/Desktop/MPCORB.DAT')
for i in range(41):
	line = infile.readline()
	#print "---",i,"---",line

def toScale( x, xmin, xmax, nsteps ):
	if (x>=xmax or x<xmin):
		return -1
	else:
		return nsteps*(x-xmin)/(xmax-xmin)

data=zeros((3,maxlines))

nsteps = 100
max_e = pi/8
max_i = pi/8

datahist = zeros((nsteps,nsteps))

print shape(datahist) 

for i in range(maxlines):
	line = infile.readline()
	#print line
	pl = planet_mpcorb(line)
	data[0,i]=pl.radius
	ix = toScale(pl.orbital_elements[1],0,max_e,nsteps )
	if (ix>0):
		datahist[ix]+=1; 
infile.close()

plot( sz );


show()
