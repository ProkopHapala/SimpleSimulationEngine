#!/usr/bin/python

import numpy as np
import sys

fname_in  = sys.argv[1]
fname_out = sys.argv[1]

def toZeroMinMax( bas ):
	xmin=bas[:,1].min(); xmax=bas[:,1].max(); bas[:,1]-=(xmin+xmax)*0.5
	ymin=bas[:,2].min(); ymax=bas[:,2].max(); bas[:,2]-=(ymin+ymax)*0.5
	zmin=bas[:,3].min(); zmax=bas[:,3].max(); bas[:,3]-=(zmin+zmax)*0.5

def toZeroAverage( bas ):
	invN = 1.0/len(bas)
	xsum=bas[:,1].sum(); bas[:,1]-=xsum*invN
	ysum=bas[:,2].sum(); bas[:,2]-=ysum*invN
	zsum=bas[:,3].sum(); bas[:,3]-=zsum*invN

# ======= main

bas=np.genfromtxt(fname_in, skip_header=1)
#toZeroMinMax( bas )
toZeroAverage( bas )
np.savetxt( fname_out,bas, header="%i"%len(bas), comments='', fmt='%i %3.6f   %3.6f   %3.6f' )


