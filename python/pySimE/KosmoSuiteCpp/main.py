import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_void_p
import ctypes
import os

from .. import utils


'''
name='KosmoSuite'

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../../'+'/cpp/libs/KosmoSuiteCpp')
print "KosmoSuiteCpp LIB_PATH     = ", LIB_PATH
print "KosmoSuiteCpp LIB_PATH_CPP = ", LIB_PATH_CPP 

if utils.recompile: 
	utils.compile_lib( name, path = LIB_PATH_CPP )

lib    = ctypes.CDLL( LIB_PATH_CPP+"/lib"+name+ utils.ext )
'''



def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../../'+'/cpp/Build/libs/KosmoSuite')
recompile(LIB_PATH_CPP)
lib = ctypes.CDLL( LIB_PATH_CPP+"/libKosmoSuite.so" )

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')


# ========= Nbody

# void set_Nbody( int n, double * mass, double * poss, double * vs ){
lib.nbody_setup.argtypes   = [c_int, array1d, array2d, array2d, array2d ]
lib.nbody_setup.restype    = None
def nbody_setup( mass, poss, vs, errs ):
	lib.nbody_setup( len(mass), mass, poss, vs, errs )

# void Nbody( double dtstart, int nstep, double * tsIn, double * tsOut, double * poss, double * vs  ){
lib.nbody_run.argtypes   = [c_double, c_double, c_double, c_int, array1d, array1d, array3d, array3d ]
lib.nbody_run.restype    = None
def nbody_run( dt_start, dt_min, dt_max, tsIn, tsOut, poss, vs ):
	lib.nbody_run( dt_start, dt_min, dt_max, len(tsIn), tsIn, tsOut, poss, vs )


# ========= ShipAccel

# void shipAccel_setup( double GM, double accel, double ivErr, double irErr, double * poss, double * vs ){
lib.shipAccel_setup.argtypes   = [ c_double, c_double, c_double, c_double, array1d, array1d,  ]
lib.shipAccel_setup.restype    = None
def shipAccel_setup( GM, accel, rErr, vErr, poss, vs ):
	lib.shipAccel_setup( GM, accel, 1.0/rErr, 1.0/vErr, poss, vs )

# shipAccel_run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
lib.shipAccel_run.argtypes   = [c_double, c_double, c_double, c_int, array1d, array1d, array2d, array2d ]
lib.shipAccel_run.restype    = None
def shipAccel_run( dt_start, dt_min, dt_max, tsIn, tsOut, poss, vs ):
	lib.shipAccel_run( dt_start, dt_min, dt_max, len(tsIn), tsIn, tsOut, poss, vs )


# ========= ElMag

# void elmag_sample ( int m, Vec3d * wheres, Vec3d * Bs, int n, Vec3d * ps, const Vec3d& dIs ){
lib.elmag_sample.argtypes   = [ c_int, array2d, array2d, c_int, array2d, array2d ]
lib.elmag_sample.restype    = None
def elmag_sample( wheres, Bs, ps, dIs ):
	lib.elmag_sample( len(wheres), wheres, Bs, len(ps), ps, dIs )

# void elmag_setup( double charge, double mass,  double ivErr, double irErr, double * poss, double * vs ){
lib.elmag_setup.argtypes   = [ c_double, c_double, c_double, c_double, array1d, array1d, c_int, array2d, array2d ]
lib.elmag_setup.restype    = None
def elmag_setup( charge, mass, rErr, vErr, poss, vs, ps, dIs ):
	lib.elmag_setup( charge, mass, 1.0/rErr, 1.0/vErr, poss, vs,  len(ps), ps, dIs   )

# void elmag_run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
lib.elmag_run.argtypes   = [c_double, c_double, c_double, c_int, array1d, array1d, array2d, array2d ]
lib.elmag_run.restype    = None
def elmag_run( dt_start, dt_min, dt_max, tsIn, tsOut, poss, vs ):
	lib.elmag_run( dt_start, dt_min, dt_max, len(tsIn), tsIn, tsOut, poss, vs )


# ========= FissionPulse

# void set_params( double mass, double alphaDOF, double total_crossection, double fission_share, double back_scatter_share, double generation_rate, double spontaneous_rate, double E_fission, double nmult, double nmult_spontal){
lib.FissionPulse_set_params.argtypes   = [ c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double ]
lib.FissionPulse_set_params.restype    = None
def FissionPulse_set_params( mass=30.0, alphaDOF=1.5, total_crossection=5.00*1e-28, fission_share=0.1, back_scatter_share=0.4, generation_rate=1/1e-8, spontaneous_rate=0.239*22.0/6.02214085774e+23, E_fission=3.20435324e-11, nmult=3.1231, nmult_spontal=2.9 ):
	lib.FissionPulse_set_params( mass, alphaDOF, total_crossection, fission_share, back_scatter_share, generation_rate, spontaneous_rate, E_fission, nmult, nmult_spontal )


# void set_initial( double R0, double v0, double Nf0, double Q0 ){
lib.FissionPulse_set_initial.argtypes   = [ c_double, c_double, c_double, c_double, c_double ]
lib.FissionPulse_set_initial.restype    = None
def FissionPulse_set_initial( R0=0.1, v0=-2e+4, Nf0=(5.5/0.239)*6.02214085774e+23, Nn0=0.0, Q0=0.0  ):
	lib.FissionPulse_set_initial( R0, v0, Nf0, Nn0, Q0  )

#void run_fixStep( double dt, int nstep, double * buff ){
lib.FissionPulse_run_fixStep.argtypes   = [ c_double, c_int,  array2d ]
lib.FissionPulse_run_fixStep.restype   = None
def FissionPulse_run_fixStep( buff, dt = 0.5e-9 ):
	lib.FissionPulse_run_fixStep( dt, len(buff), buff )


# void set_initial( double time_trigger, double R_trigger, double Nn_spark ){
lib.FissionPulse_set_trigger.argtypes   = [ c_double, c_double, c_double ]
lib.FissionPulse_set_trigger.restype    = None
def FissionPulse_set_trigger( time_trigger=1e+300, R_trigger=0.05, Nn_spark=1e+18 ):
	lib.FissionPulse_set_trigger( time_trigger, R_trigger, Nn_spark )


# ========= SpaceLaunchODE
# simple example is in /home/prokop/Dropbox/MyDevSW/ctypes/PassStruct
# this things should be automatized
#see https://github.com/davidjamesca/ctypesgen/tree/master/ctypesgencore
#    http://cffi.readthedocs.io/en/latest/
#    http://stackoverflow.com/questions/14534998/tool-to-convert-c-structure-to-ctypes-structure

class Vec3d(ctypes.Structure):
     _fields_ = [
     ("x",    c_double),
     ("y",    c_double),
     ("z",    c_double)
     ]

class Launch(ctypes.Structure):
     _fields_ = [
     ("n",         c_int),
     ("thetaCPs ", c_void_p),
     ("thrustCPs", c_void_p),
     ("dirCPs",    c_void_p),
     ("uT",        c_double),
     ("inv_uT",    c_double),
     ("tMax",      c_double),
     ("vTarget",   c_double),
     ("hTarget",   c_double)
     ]

class Aerodynamics(ctypes.Structure):
     _fields_ = [
     ("n",       c_int),
     ("dvM",     c_double),
     ("inv_dvM", c_double),
     ("vMax",    c_double),
     ("Cd_CPs",  c_void_p)
     ]

class Atmosphere(ctypes.Structure):
     _fields_ = [
     ("n",        c_int),
     ("dh",       c_double),
     ("inv_dh",   c_double),
     ("hmax",     c_double),
     ("rho_CPs",  c_void_p),
     ("rho0",     c_double),
     ("zrate",    c_double)
     ]

class Rocket(ctypes.Structure):
     _fields_ = [
     ("vexhaust",    c_double),
     ("dm_F",        c_double),
     ("mass_initial",c_double),
     ("mass_empty",  c_double),
     ("thrust_full", c_double),
     ("AeroArea",    c_double),
     ]

class Planet(ctypes.Structure):
     _fields_ = [
     ("R",   c_double),
     ("GM",  c_double),
     ("pos", Vec3d)
     ]

class LogTrig(ctypes.Structure):
     _fields_ = [
     ("on",      c_bool),
     ("dt_trig", c_double),
     ("t_trig",  c_double),
     ("i",       c_int),
     ]

# void SpaceLaunchODE_init( Planet *planet_, Rocket *rocket_, Aerodynamics *aero_, Atmosphere *atmosphere_){
lib.SpaceLaunchODE_init.argtypes   = [ Planet, Rocket, Aerodynamics, Atmosphere ]
lib.SpaceLaunchODE_init.restype    = None
#def SpaceLaunchODE_init( planet, rocket, aero, atmosphere ):
#	lib.SpaceLaunchODE_init( planet, rocket, aero, atmosphere )

# int SpaceLaunchODE_run( int nLogMax, int nMaxIters, Launch *launch_, LogTrig *logTrig_, double * outbuff ){
lib.SpaceLaunchODE_run.argtypes   = [ c_int, c_int, Launch, LogTrig, ctypes.c_void_p ]
lib.SpaceLaunchODE_run.restype    = c_int
#def SpaceLaunchODE_run( ):
#	lib.SpaceLaunchODE_run(  )




