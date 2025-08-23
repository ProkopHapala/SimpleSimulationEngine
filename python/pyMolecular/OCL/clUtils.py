import sys
import os
import numpy as np
import ctypes
import pyopencl as cl
# import pyopencl.array as cl_array
# import pyopencl.cltypes as cltypes
# import matplotlib.pyplot as plt
import time

bytePerFloat = 4

FFT = None


import pyopencl as cl

def get_nvidia_device( what="nvidia"):
    platforms = cl.get_platforms()
    for platform in platforms:
        devices = platform.get_devices()
        for device in devices:
            if what in device.name.lower():
                # Create the OpenCL context and command queue
                ctx = cl.Context([device])
                queue = cl.CommandQueue(ctx)
                # Print information about the selected device
                print(f"Selected device: {device.name}")
                get_cl_info(device)
                return ctx, queue
    # If no NVIDIA device is found, return None
    print("pyOepnCL error: No {what} device found.")
    return None, None


def try_load_clFFT():
    global FFT
    if FFT is None:
        from gpyfft.fft import FFT as FFT_
        FFT = FFT_

def make_inds_pbc(n):
    return np.array([
    [ 0, 1,   2,   3   ],
    [ 0, 1,   2,   3-n ],
    [ 0, 1,   2-n, 3-n ],
    [ 0, 1-n, 2-n, 3-n ]], dtype=np.int32 );

def roundup_global_size(global_size, local_size):
    remainder = global_size % local_size
    if remainder == 0: return global_size
    return global_size + (local_size - remainder)

def roundup_global_size_3d( global_size,  local_size):
    return (
         roundup_global_size(global_size[0], local_size[0]),
         roundup_global_size(global_size[1], local_size[1]),
         roundup_global_size(global_size[2], local_size[2])
     )

def get_cl_info( device ):
    '''
    also use:
        nvidia-smi -q
    '''
    # Example for NVIDIA GPUs:
    # Pascal: 128 CUDA cores per SM
    # Turing: 64 CUDA cores per SM
    # Ampere: 128 CUDA cores per SM  10496 cores / 82 compute units

    print(f"Device Name: {device.name}")
    print(f"Max Compute Units: {device.max_compute_units}")
    print(f"Max Work Group Size: {device.max_work_group_size}")
    print(f"Global Memory Size: {device.global_mem_size / (1024*1024)} MB")
    print(f"Local Memory Size: {device.local_mem_size / 1024} KB")
    print(f"Max Clock Frequency: {device.max_clock_frequency} MHz")
              
    # Get local memory characteristics - handle missing characterize module gracefully
    try:
        granularity      = cl.characterize.local_memory_access_granularity(device)
        bank_count       = cl.characterize.local_memory_bank_count(device)
        usable_local_mem = cl.characterize.usable_local_mem_size(device)
        # Print results
        print(f"Local Memory Access Granularity: {granularity} bytes")
        print(f"Number of Local Memory Banks: {bank_count}")
        print(f"Usable Local Memory Size: {usable_local_mem} bytes")
    except AttributeError as e:
        print(f"Note: PyOpenCL characterize module not available. Some device info will not be displayed.")

    # Retrieve various characteristics
    try:
        fast_math_options = cl.characterize.get_fast_inaccurate_build_options(device)
        simd_group_size   = cl.characterize.get_simd_group_size(device, 4)  # Assuming float4
        double_support    = cl.characterize.has_amd_double_support(device)
        #src_cache_support = cl.characterize.has_src_build_cache(device)

        # Print additional information
        print(f"SIMD Group Size: {simd_group_size}")
        print(f"Has AMD Double Support: {double_support}")
        print(f"Fast Math Options: {fast_math_options}")
    except AttributeError as e:
        print(f"Note: Additional PyOpenCL characterization info not available.")

def local_memory_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return device.local_mem_size /( sp_per_cu/local_size )

def local_memory_float_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/4

def local_memory_float4_per_workgroup( device, local_size=32, sp_per_cu=128 ):
    return local_memory_per_workgroup( device, local_size=local_size, sp_per_cu=sp_per_cu )/(4*4)

class GridShape:
    def __init__(self, ns=None, dg=(0.0,0.0,0.0), lvec=None, Ls=None, g0=(0.0,0.0,0.0) ):
        if lvec is None: lvec = np.array( [ [Ls[0],0.,0.], [0.,Ls[1],0.], [0.,0.,Ls[2]] ] )
        if Ls   is None: Ls = (lvec[0][0],lvec[1][1],lvec[2][2])
        if ns   is None: ns = (int((Ls[0]/dg[0])+0.5),int((Ls[1]/dg[1])+0.5),int((Ls[2]/dg[2])+0.5))
        self.ns   = ns
        self.nxyz = np.prod(ns)
        self.Ls   = Ls
        self.g0 = g0
        self.dg   = dg
        self.lvec = lvec
        self.dV   = dg[0]*dg[1]*dg[2] 
        self.V    = self.dV*self.nxyz
        

    def __str__(self):
        return f"GridShape(ns={self.ns}, Ls={self.Ls}, pos0={self.pos0}, dg={self.dg}, lvec={self.lvec})"

class GridCL:
    def __init__(self, gsh : GridShape ):
        self.nxyz = np.int32( gsh.ns[0]*gsh.ns[1]*gsh.ns[2] )
        self.ns = np.array( gsh.ns+(self.nxyz,), dtype=np.int32   )
        self.g0 = np.array( gsh.g0+(0.,),   dtype=np.float32 )
        self.dg = np.array( gsh.dg+(0.,),   dtype=np.float32 )
        self.a  = np.array( ( *gsh.lvec[0], 0.), dtype=np.float32 )
        self.b  = np.array( ( *gsh.lvec[1], 0.), dtype=np.float32 )
        self.c  = np.array( ( *gsh.lvec[2], 0.), dtype=np.float32 )
        self.da = np.array( ( *(gsh.lvec[0]/gsh.ns[0]), 0.), dtype=np.float32 )
        self.db = np.array( ( *(gsh.lvec[1]/gsh.ns[1]), 0.), dtype=np.float32 )
        self.dc = np.array( ( *(gsh.lvec[2]/gsh.ns[2]), 0.), dtype=np.float32 )
        
    def __str__(self):
        return f"GridShape(ng={self.ng}, g0={self.g0}, dg={self.dg}\n   La={self.La},Lb={self.Lb},Lc={self.Lc}\n    da={self.da},db={self.db},dc={self.dc})"
