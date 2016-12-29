#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def getSomeOrtho( axis ):
    axis  = axis/np.linalg.norm(axis)
    np.dot([])
    ys
    return np.array([left,up,fw])

def makeRotMat( fw, up ):
    fw   = fw/np.linalg.norm(fw)
    up   = up - fw*np.dot(up,fw)
    up   = up/np.linalg.norm(up)
    left = np.cross(fw,up)
    left = left/np.linalg.norm(left) 
    return np.array([left,up,fw])
    
def plotLines( p1s, p2s, clr=None, cam ):
    return None

def plotCurve( ps, cam ):
    ps_ = np.dot( cam, ps )
    plt.plot( ps[0], ps[1] )

def plotSpots( sz, cam ):
    return None
    
    
def circle( p0, R, axis, n=32 ):    
    ts    = np.linspace(10,10,R)
    
    ca = R * np.sin(ts)
    sa = R * np.sin(ts)
    return 
    
    
    
    
    
    
    
if __name__ == "__main__":
    
    np.sin(ts)
