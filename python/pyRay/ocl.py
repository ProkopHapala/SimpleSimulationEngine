#!/usr/bin/python

import os
import pyopencl as cl
import numpy    as np 

from common import *

plats      = cl.get_platforms()
ctx        = cl.Context(properties=[(cl.context_properties.PLATFORM, plats[0])], devices=None)       
queue      = cl.CommandQueue(ctx)
cl_program = None

'''
def loadProgram(fname, ctx=ctx, queue=queue):
    f       = open(fname, 'r')
    fstr    = "".join(f.readlines())
    program = cl.Program(ctx, fstr).build()
    return program
'''

def make_program(src):
    global cl_program
    cl_program = cl.Program(ctx, src ).build()
    return cl_program

def getRays( cam, fw0=-4.0, t0=8.0, tg=(0.25,0.25), nsz=screen_resolution ):
    #print cam
    xs     = np.linspace(-tg[1],tg[1],nsz[1])
    ys     = np.linspace(-tg[0],tg[0],nsz[0])
    Xs,Ys  = np.meshgrid(xs,ys)
    #print cam.shape
    #print nsz,Xs.shape,Ys.shape
    rd     = np.zeros(nsz+[4,], dtype=np.float32 )
    ro     = np.zeros(nsz+[4,], dtype=np.float32 )
    rd[:,:,0] = Xs*cam[0][0] + Ys*cam[1][0] + cam[2][0]
    rd[:,:,1] = Xs*cam[0][1] + Ys*cam[1][1] + cam[2][1]
    rd[:,:,2] = Xs*cam[0][2] + Ys*cam[1][2] + cam[2][2]
    rd[:,:,3] = np.sqrt( rd[:,:,0]**2 + rd[:,:,1]**2 + rd[:,:,2]**2 )
    rd[:,:,0:3] /= rd[:,:,3][:,:,None]
    ro[:,:,0:3]   = cam[None,None,2]*fw0
    #ro[:,:,0:3] += rd[:,:,0:3]*t0*t0
    #print ro
    return rd, ro

def prep_rayTraceBasic( rd, ro ):
    #print "rd = ", rd
    #print "ro = ", ro
    mf      = cl.mem_flags
    cl_rd   = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=rd )
    cl_ro   = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=ro )
    cl_hitp = cl.Buffer(ctx, mf.WRITE_ONLY, ro.nbytes )
    cl_hitn = cl.Buffer(ctx, mf.WRITE_ONLY, ro.nbytes )
    kargs = ( cl_rd, cl_ro, cl_hitp, cl_hitn )
    return kargs 
    	
def run_rayTraceBasic( kargs, nsz=screen_resolution, local_size=(16,), max_depth=20.0 ):
    #print "run opencl kernel ..."
    global_size = (nsz[0]*nsz[1],)
    hitp        = np.zeros( nsz+[4,] , dtype=np.float32 )
    hitn        = np.zeros( nsz+[4,] , dtype=np.float32 )
    args_       = (np.float32(max_depth),)
    cl_program.rayTrace_basic( queue, global_size, local_size, *(kargs+args_) )
    cl.enqueue_copy          ( queue, hitp, kargs[2] );
    cl.enqueue_copy          ( queue, hitn, kargs[3] );
    queue.finish()
    #print "... opencl kernel DONE"
    return hitp,hitn
    
    
