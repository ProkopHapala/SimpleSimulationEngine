
from PyKEP import *
from pylab import *

maxlines=1000000

path='D:\\know\\Technika\\SpaceTech\\Asteroids\\Datafiles\\'
sortdir='sorted\\'
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

# =========== Sort just by a ( semimajor axis )

outfs=[None]*na
for i in range(maxlines):
	line = infile.readline()
	if (line is None) or (len(line)<10):
		break
	if ((i % 1000) ==0 ):
		print i
	pl = planet_mpcorb(line)
	a = pl.orbital_elements[0]/AU
	ia = toScale(a, min_a,max_a,na )
	if ( ia>0  ):
		if ( outfs[ia] is None ):
			fname = 'tile_a%02i.dat' %ia
			print "created file ",fname
			outfs[ia] = open( path+sortdir+fname, 'w' )
		outfs      [ia].write(line)
infile.close()
for ia in range(na):
	if ( outfs[ia] is not None ):
		outfs[ia].close()

print ' Sort just by a .... DONE ! '

# ============ Sort by a, Lx, Ly

hist_all  = zeros((na,nL,nL))

countmapf = open(path+sortdir+'count_a_Lx_Ly.dat' , 'w')
countmapf.write('    na= %5i  nLx= %5i nLy= %5i  \n' %( na, nL, nL )) 
countmapf.write(' min_a= %16.8f min_Lx= %16.8f min_Ly= %16.8f  \n' %( min_a, -max_L, -max_L )) 
countmapf.write(' max_a= %16.8f max_Lx= %16.8f max_Ly= %16.8f  \n' %( max_a,  max_L,  max_L )) 
for ia in range(na):
	outfs = []
	for ix in range(nL):
		outfs.append( [None]*nL)
	inname = 'tile_a%02i.dat' %ia
	print " opening ... ",inname
	try:
		infile  = open(path+sortdir+inname,'r')
		for i in range(maxlines):
			line = infile.readline()
			if (line is None) or (len(line)<10):
				break
			pl = planet_mpcorb(line)
			si = sin(pl.orbital_elements[2])
			cW = cos(pl.orbital_elements[3])
			sW = sin(pl.orbital_elements[3])
			x  = si*cW
			y  = si*sW
			ix = toScale(x,-max_L,max_L,nL )
			iy = toScale(y,-max_L,max_L,nL )
			#print ix,iy
			if ( (ix>0) and (iy>0) ):
				if ( outfs[ix][iy] is None ):
					fname = 'tile_a%02i_Lx%02i_Ly%02i.dat' %(ia,ix,iy)
					print "created file ",fname
					outfs[ix][iy] = open( path+sortdir+fname, 'w' )
					outfs[ix][iy].write( ' a= %16.8f Lx= %16.8f Ly= %16.8f   \n' %( invScale(ia,min_a,max_a,na), invScale(iy,-max_L,max_L,nL), invScale(iy,-max_L,max_L,nL)   ) )
				outfs      [ix][iy].write(line)
				hist_all   [ia,ix,iy]+=1; 
		infile.close()
		for ix in range(nL):
			for iy in range(nL):
				if outfs[ix][iy] is not None:
					outfs[ix][iy].close()
					countmapf.write(' %5i   %5i   %5i   %10i  \n' %( ia,  ix,  iy, hist_all[ia,ix,iy] )) 
	except:
		pass
countmapf.close()

print ' Sort by a, Lx, Ly .... DONE ! '

show()
