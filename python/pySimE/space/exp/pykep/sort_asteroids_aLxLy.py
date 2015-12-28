
from PyKEP import *
from pylab import *

maxlines=1000000

path='D:\\know\\Technika\\SpaceTech\\Asteroids\\Datafiles\\'
sordir='sorted\\'
#path='/home/prokop/Desktop/'
infile  = open(path+'MPCORB.DAT','r')


for i in range(41):
	line = infile.readline()
	#print "---",i,"---",line

def toScale( x, xmin, xmax, nsteps ):
	if (x>=xmax or x<xmin):
		return -1
	else:
		return int(nsteps*(x-xmin)/(xmax-xmin))
		
def invScale( ix, xmin, xmax, nsteps ):
		return (float(ix)*(xmax-xmin)/nsteps)+xmin



na = 12; min_a =2.1; max_a=3.3;
nL = 10; max_L  = 0.3

hist_all  = zeros((na,nL,nL))

outfs=[]
for ia in range(na):
	atemp = []
	for ix in range(nL):
		atemp.append( [None]*nL)
	outfs.append(atemp)

for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	pl = planet_mpcorb(line)
	a = pl.orbital_elements[0]/AU
	si = sin(pl.orbital_elements[2])
	cW = cos(pl.orbital_elements[3])
	sW = sin(pl.orbital_elements[3])
	x  = si*cW
	y  = si*sW
	ix = toScale(x,-max_L,max_L,nL )
	iy = toScale(y,-max_L,max_L,nL )
	ia = toScale(a, min_a,max_a,na )
	#print ix,iy
	if ( (ia>0) and (ix>0) and (iy>0) ):
		if ( outfs[ia][ix][iy]==None ):
			fname = 'tile_a%02i_Lx%02i_Ly%02i.dat' %(ia,ix,iy)
			print "created file ",fname
			outfs[ia][ix][iy] = open( path+sordir+fname, 'w' )
			outfs[ia][ix][iy].write( ' a= %16.8f Lx= %16.8f Ly= %16.8f   \n' %( invScale(ia,min_a,max_a,na), invScale(iy,-max_L,max_L,nL), invScale(iy,-max_L,max_L,nL)   ) )
		outfs      [ia][ix][iy].write(line)
		hist_all   [ia,ix,iy]+=1; 
infile.close()

countmapf = open(path+'count_a_Lx_Ly.dat' , 'w')
countmapf.write('    na= %5i  nLx= %5i nLy= %5i  \n' %( na, nL, nL )) 
countmapf.write(' min_a= %16.8f min_Lx= %16.8f min_Ly= %16.8f  \n' %( min_a, -max_L, -max_L )) 
countmapf.write(' max_a= %16.8f max_Lx= %16.8f max_Ly= %16.8f  \n' %( max_a,  max_L,  max_L )) 
for ia in range(na):
	atemp = []
	for ix in range(nL):
		for iy in range(nL):
			if outfs[ia][ix][iy]!=None:
				outfs[ia][ix][iy].close()
				countmapf.write(' %5i   %5i   %5i   %10i  \n' %( ia,  ix,  iy, hist_all[ia,ix,iy] )) 
countmapf.close()

print "done"

show()
