#!/usr/bin/python

import os
import numpy as np


defines={
"MAX_STEPS": 256, 
"HIT_PREC" : 0.0001, 
}

PACKAGE_PATH = os.path.dirname( os.path.realpath( __file__ ) ); print PACKAGE_PATH
CL_PATH      = os.path.normpath( PACKAGE_PATH + '/cl' )

default_up        = np.array([0.0,1.0,0.0])
screen_resolution = [256,256]

def getCamMat( fw, up=default_up ):
    fw/= np.sqrt( np.dot(fw,fw) )
    up = up -  fw*np.dot(up,fw)
    up/= np.sqrt( np.dot(up,up) )
    lf = np.cross(up,fw)
    return np.array([lf,up,fw])



