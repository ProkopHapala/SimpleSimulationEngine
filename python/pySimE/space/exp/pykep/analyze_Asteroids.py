
from PyKEP import *
from pylab import *

maxlines=10000

iline = 0
infile  = open('/home/prokop/Desktop/MPCORB.DAT')
for i in range(41):
	line = infile.readline()
	#print "---",i,"---",line



data=zeros((3,maxlines))
for i in range(maxlines):
	line = infile.readline()
	#print line
	pl = planet_mpcorb(line)
	#print pl.orbital_elements
	#print pl.radius, log(pl.radius)
	data[0,i]=pl.radius
	data[1,i]=pl.orbital_elements[1]
	data[2,i]=pl.orbital_elements[2]	
infile.close()


#scatter(data[1], data[2], s=log(data[0]) )
scatter(data[1], data[2], s=1 )
show()
