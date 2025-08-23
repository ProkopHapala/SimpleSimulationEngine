import sys
import os
import numpy as np
import re

import pyopencl as cl
# from . import clUtils as clu
import pyopencl.array as cl_array
import pyopencl.cltypes as cltypes
import matplotlib.pyplot as plt
import time

from . import clUtils as clu
from .MMFF import MMFF
from .OpenCLBase import OpenCLBase

REQ_DEFAULT = np.array([1.7, 0.1, 0.0, 0.0], dtype=np.float32)  # R, E, Q, padding

verbose=False

def pack(iSys, source_array, target_buffer, queue):
    offset = iSys * source_array.size * source_array.dtype.itemsize
    cl.enqueue_copy(queue, target_buffer, source_array, offset=offset)

def copy(source, target, queue, iSys):
    pack(iSys, source, target, queue)

def copy_add(source, source_add, target, offset, queue):
    pass

def mat3_to_cl(mat3_np):
    return mat3_np.flatten().astype(np.float32)

def vec3_to_cl(vec3_np):
    return np.append(vec3_np, 0.0).astype(np.float32)

class MolecularDynamics(OpenCLBase):
    """
    Class for molecular dynamics simulations using OpenCL.
    
    This class inherits from OpenCLBase and implements specific functionality
    for molecular dynamics simulations using the relax_multi_mini.cl kernel.
    """
    
    def __init__(self, nloc=32, perBatch=10):
        # Initialize the base class
        super().__init__(nloc=nloc, device_index=0)
        
        # Load the OpenCL program
        base_path = os.path.dirname(os.path.abspath(__file__))
        rel_path = "cl/relax_multi_mini.cl"
        if not self.load_program(rel_path=rel_path, base_path=base_path, bPrint=False):
            exit(1)
        
        # Initialize other attributes that will be set in realloc
        self.nSystems       = 0
        self.mmff_list      = []
        self.MD_event_batch = None
        self.perBatch       = perBatch
        self.nstep          = 1

    def realloc(self, mmff, nSystems=1 ):
        """
        Reallocate buffers for the given number of systems based on the MMFF template.
        """
        # Store dimensions explicitly to avoid reference issues
        print(f"MolecularDynamics::realloc() natoms={mmff.natoms}, nvecs={mmff.nvecs}, nnode={mmff.nnode}")
        self.nSystems = nSystems
        self.mmff_list = [mmff] * nSystems  # Assuming all systems use the same MMFF parameters
        self.allocate_cl_buffers(mmff)
        self.allocate_host_buffers()

    def allocate_host_buffers(self):
        self.atoms  = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)
        self.aforce = np.zeros((self.nSystems, self.nvecs, 4), dtype=np.float32)

    def allocate_cl_buffers(self, mmff):
        """
        Allocates OpenCL buffers based on the MMFF template and number of systems.
        Includes all buffers required by the runMD kernel.
        """
        nSystems = self.nSystems
        natoms = mmff.natoms
        nvecs  = mmff.nvecs
        nnode  = mmff.nnode
        ncap   = mmff.ncap
        ntors  = mmff.ntors
        nbkng  = nvecs
        nPBC   = mmff.nPBC
        npbc   = mmff.npbc

        self.nDOFs = (natoms,nnode)

        self.natoms = natoms
        self.nvecs  = nvecs
        self.nnode  = nnode
        self.ncap   = ncap
        self.ntors  = ntors
        self.nbkng  = nbkng
        self.nPBC   = nPBC
        self.npbc   = npbc
        
        print(f"MolecularDynamics::allocate_cl_buffers(): nSystems: {nSystems}  natoms: {natoms}  nvecs: {nvecs} nnode: {nnode} ncap: {ncap}  ntors: {ntors}  nbkng: {nbkng}")
        
        if nSystems <= 0 or natoms <= 0 or nvecs <= 0 or nnode <= 0:
            raise ValueError(f"Invalid dimensions for buffer allocation: nSystems={nSystems}, natoms={natoms}, nvecs={nvecs}, nnode={nnode}")
        
        float_size = np.float32().itemsize
        int_size = np.int32().itemsize
        mat3_size = 3 * 4 * float_size  # 3x3 matrix
        
        mf = cl.mem_flags
        
        # Dynamical variables
        self.create_buffer('apos',     nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('aforce',   nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('avel',     nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('fneigh',   nSystems * nnode * 4 * 2 * float_size, mf.READ_WRITE)
        self.create_buffer('cvf',      nSystems * nvecs * 4 * float_size, mf.READ_WRITE)
        
        # Neighbor lists
        self.create_buffer('neighs',    nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('neighCell', nSystems * natoms * 4 * int_size, mf.READ_ONLY)
        self.create_buffer('bkNeighs',  nSystems * nbkng * 4 * int_size, mf.READ_ONLY)
        
        # Force field parameters
        self.create_buffer('REQs',     nSystems * natoms * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('apars',    nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('bLs',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('bKs',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('Ksp',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('Kpp',      nSystems * nnode * 4 * float_size, mf.READ_ONLY)
        
        # System parameters
        self.create_buffer('lvecs',    nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('ilvecs',   nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('pbc_shifts', nSystems * npbc * 4 * float_size, mf.READ_ONLY)
        
        # MD parameters and constraints
        self.create_buffer('constr',   nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('constrK',  nSystems * natoms * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('MDparams', nSystems * 4 * float_size, mf.READ_ONLY)
        self.create_buffer('TDrives',  nSystems * 4 * float_size, mf.READ_ONLY)
        
        # System interactions
        self.create_buffer('bboxes',   nSystems * mat3_size, mf.READ_ONLY)
        self.create_buffer('sysneighs',nSystems * int_size, mf.READ_ONLY)
        self.create_buffer('sysbonds', nSystems * 4 * float_size, mf.READ_ONLY)
        
        # Grid force field parameters (scalar)
        # self.kernel_params['GFFParams'] = np.zeros(4, dtype=np.float32)
        # self.kernel_params['bSubtractVdW'] = np.int32(0)

    def pack_system(self, iSys, mmff):
        """Packs data from an MMFF instance into GPU buffers for a specific system index."""
        #print("pack_system() iSys=%d" % iSys)
        nvecs   = mmff.nvecs
        natoms  = mmff.natoms
        nnode   = mmff.nnode
        float4_size = 4 * np.float32().itemsize
        int4_size   = 4 * np.int32().itemsize

        offset_atoms = iSys * nvecs * float4_size
        self.toGPU('apos',      mmff.apos.astype(np.float32).flatten(), byte_offset=offset_atoms)
        self.toGPU('aforce',    mmff.fapos.astype(np.float32).flatten(), byte_offset=offset_atoms)

        offset_REQs = iSys * natoms * float4_size
        self.toGPU('REQs',      mmff.REQs.astype(np.float32).flatten(), byte_offset=offset_REQs)
        
        offset_neighs = iSys * natoms * int4_size
        self.toGPU('neighs',    mmff.neighs.astype(np.int32).flatten(), byte_offset=offset_neighs)
        self.toGPU('neighCell', mmff.neighCell.astype(np.int32).flatten(), byte_offset=offset_neighs)
        
        offset_apars = iSys * nnode * float4_size
        self.toGPU('apars',     mmff.apars.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bLs',       mmff.bLs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('bKs',       mmff.bKs.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Ksp',       mmff.Ksp.astype(np.float32).flatten(), byte_offset=offset_apars)
        self.toGPU('Kpp',       mmff.Kpp.astype(np.float32).flatten(), byte_offset=offset_apars)
        #print("pack_system() iSys=%d" % iSys, "offset_atoms=%d" % offset_atoms, "offset_REQs=%d" % offset_REQs, "offset_neighs=%d" % offset_neighs, "offset_apars=%d" % offset_apars)
        #print("pack_system() iSys=%d" % iSys, "atoms.nbytes=%d" % mmff.apos.nbytes, "REQs.nbytes=%d" % mmff.REQs.nbytes, "neighs.nbytes=%d" % mmff.neighs.nbytes, "apars.nbytes=%d" % mmff.apars.nbytes)
        self.toGPU('MDparams',  np.array([mmff.dt, mmff.damp, mmff.Flimit], dtype=np.float32), byte_offset=iSys*float4_size)


    def init_with_atoms(self, na=None, atoms=None, REQs=None, REQ_default=REQ_DEFAULT):
        """
        Initialize MolecularDynamics directly with atom positions and REQ parameters,
        without using MMFF object.
        
        Args:
            atoms: numpy array of shape (na, 3) with atom positions
            REQs: numpy array of shape (na, 4) with REQ parameters or None for defaults
            nSystems: number of systems
            nloc: local workgroup size
            
        Returns:
            Initialized MolecularDynamics instance
        """
        float_size = np.float32().itemsize

        if na is None:
            na = len(atoms)
        mf = cl.mem_flags

        # Initialize necessary attributes
        self.nSystems = 1
        self.natoms   = na
        self.nvecs    = na   # We don't have nodes or pi-orbitals, just atoms
        self.nnode    = 0         # No node atoms
        self.nDOFs    = (na, 0)  # (natoms, nnode)
        
        # Additional required attributes for initGridFF
        self.nPBC  = (0, 0, 0)  # No periodic boundary conditions
        self.npbc  = 0
        self.nbkng = na
        self.ncap  = 0
        self.ntors = 0
        # Convert atoms to float32 if needed
        
        if atoms is None:
            self.atoms = np.zeros((na, 4), dtype=np.float32)
        else:
            self.atoms = np.asarray(atoms, dtype=np.float32)
        if REQs is None:
            self.REQs    = np.zeros((na, 4), dtype=np.float32)
            self.REQs[:, 0] = REQ_default[0]
            self.REQs[:, 1] = REQ_default[1]
            self.REQs[:, 2] = REQ_default[2]
            self.REQs[:, 3] = REQ_default[3]
        else:
            self.REQs    = np.asarray(REQs, dtype=np.float32)
        self.aforce = np.zeros((na, 4), dtype=np.float32)

        self.create_buffer('apos',   na * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('aforce', na * 4 * float_size, mf.READ_WRITE)
        self.create_buffer('REQs',   na * 4 * float_size, mf.READ_ONLY)

        self.toGPU('apos',   self.atoms  )
        self.toGPU('aforce', self.aforce )
        self.toGPU('REQs',   self.REQs   )
        self.queue.finish()

    def setup_kernels(self):
        """
        Prepares the kernel arguments for all kernels by parsing their headers.
        Also sets up the work sizes for each kernel.
        """
        print("MolecularDynamics::setup_kernels()")
        # Get all work sizes at once
        self.get_work_sizes()
        self.init_kernel_params()
                
        # Generate kernel arguments
        self.kernel_args_getMMFFf4         = self.generate_kernel_args("getMMFFf4")
        self.kernel_args_getNonBond        = self.generate_kernel_args("getNonBond")
        # --- NOTE: grid-kernels are intialized in initGridFF()
        #self.kernel_args_getNonBond_GridFF_Bspline = self.generate_kernel_args("getNonBond_GridFF_Bspline")
        #self.kernel_args_getNonBond_GridFF_Bspline_tex = self.generate_kernel_args("getNonBond_GridFF_Bspline_tex")
        self.kernel_args_updateAtomsMMFFf4 = self.generate_kernel_args("updateAtomsMMFFf4")
        self.kernel_args_cleanForceMMFFf4  = self.generate_kernel_args("cleanForceMMFFf4")
        self.kernel_args_runMD             = self.generate_kernel_args("runMD")

    def init_kernel_params(self):
        """
        Initialize a dictionary of standard kernel parameters.
        This provides default values for common parameters used in kernels.
        """
        print("MolecularDynamics::init_kernel_params()")
        # Create a dictionary to store kernel parameters
        self.kernel_params = {
            # Common dimension parameters
            'nDOFs':        np.array([self.natoms, self.nnode, 0, self.perBatch], dtype=np.int32),
            'mask':         np.array([1, 1, 1, 1],         dtype=np.int32),
            'nPBC':         np.array(self.nPBC+(0,),       dtype=np.int32),
            'GFFParams':    np.array([0.0, 0.0,  0.0, 0.0], dtype=np.float32),
            'MDparams':     np.array([0.1, 0.05, 0.0, 0.0], dtype=np.float32),
            # Common scalar parameters
            'npbc':         np.int32(self.npbc),
            'bSubtractVdW': np.int32(0),
            'grid_ns':      np.array([0,0,0,0], dtype=np.int32),
            'grid_invStep': np.array([0.0,0.0,0.0,0.0], dtype=np.float32),
            'grid_p0':      np.array([0.0,0.0,0.0,0.0], dtype=np.float32),
        }
        
    def get_work_sizes(self):
        """
        Generate standard work sizes based on current system dimensions.
        
        Returns:
            dict: Dictionary containing sz_na, sz_node, sz_nvec, sz_loc
        """
        # Default to the standard work size parameters
        self.sz_loc = (self.nloc, 1)
        self.sz_na   = (clu.roundup_global_size(self.natoms, self.nloc), self.nSystems)
        self.sz_node = (clu.roundup_global_size(self.nnode,  self.nloc), self.nSystems)
        self.sz_nvec = (clu.roundup_global_size(self.nvecs,  self.nloc), self.nSystems)
        
        # Return all sizes, let the caller decide which to use
        return {
            'sz_na':   self.sz_na,
            'sz_node': self.sz_node,
            'sz_nvec': self.sz_nvec,
            'sz_loc':  self.sz_loc
        }

    def upload_all_systems(self):
        """Uploads data for all systems to the GPU."""
        for sys_idx in range(self.nSystems):
            self.pack_system(sys_idx, self.mmff_list[sys_idx])
        #print("MolecularDynamics::upload_all_systems() DONE")

    def run_getNonBond(self):
        self.prg.getNonBond(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond)
        self.queue.finish()
    
    def run_getNonBond_GridFF_Bspline(self):
        self.prg.getNonBond_GridFF_Bspline(self.queue, self.sz_na, self.sz_loc,  *self.kernel_args_getNonBond_GridFF_Bspline)
        self.queue.finish()

    def run_getNonBond_GridFF_Bspline_tex(self):
        self.prg.getNonBond_GridFF_Bspline_tex(self.queue, self.sz_na, self.sz_loc,  *self.kernel_args_getNonBond_GridFF_Bspline_tex)
        self.queue.finish()

    def run_getMMFFf4(self):
        self.prg.getMMFFf4(self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4)
        self.queue.finish()
    
    def run_updateAtomsMMFFf4(self):
        self.prg.updateAtomsMMFFf4(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4)
        self.queue.finish()
    
    def run_cleanForceMMFFf4(self):
        self.prg.cleanForceMMFFf4(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_cleanForceMMFFf4)
        self.queue.finish()
    
    def run_runMD(self):
        self.prg.runMD(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_runMD)
        self.queue.finish()
    
    def run_sampleGrid_tex(self, apos=None, bUseTexture=False):
        if apos is not None:
            self.toGPU('apos', apos)
        if bUseTexture:
            self.prg.sampleGrid_tex(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_sampleGrid_tex)
        else:
            self.prg.sampleGrid(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_sampleGrid)
        self.fromGPU('aforce', self.aforce)
        self.queue.finish()
        return self.aforce.reshape(-1, 4).copy()

    def run_MD_py(self, nsteps, use_gridff=False):
        """Run molecular dynamics simulation..."""
        for i in range(nsteps):
            if use_gridff and self.has_gridff:
                if self.use_texture:
                    self.prg.getNonBond_GridFF_Bspline_tex(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline_tex)
                else:
                    self.prg.getNonBond_GridFF_Bspline(self.queue, self.sz_na, self.sz_loc, *self.kernel_args_getNonBond_GridFF_Bspline)
            else:
                self.prg.getNonBond       (self.queue, self.sz_na, self.sz_loc,*self.kernel_args_getNonBond)
            self.prg.getMMFFf4        (self.queue, self.sz_node, self.sz_loc, *self.kernel_args_getMMFFf4)
            self.prg.updateAtomsMMFFf4(self.queue, self.sz_nvec, self.sz_loc, *self.kernel_args_updateAtomsMMFFf4)
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        self.queue.finish()
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)

    def download_results(self):
        self.fromGPU('apos',   self.atoms)
        self.fromGPU('aforce', self.aforce)
        return self.atoms.reshape(-1, 4), self.aforce.reshape(-1, 4)
    
    def initGridFF(self, grid_shape, bspline_data, grid_p0, grid_step, use_texture=False, r_damp=0.0, alpha_morse=0.0, bKernels=True):
        """Initialize GridFF with B-spline data"""
        
        #grid_shape = grid_shape[::-1] #.copy()

        print("MolecularDynamics::initGridFF() grid_shape: ", grid_shape)
        self.has_gridff = True
        self.use_texture = use_texture
        
        # 1. Ensure kernel_params exists
        if not hasattr(self, 'kernel_params'):
            self.init_kernel_params()
            
        # 2. Set grid parameters
        self.kernel_params.update({
            'grid_ns':      np.array([*grid_shape ,0],                   dtype=np.int32),
            'grid_invStep': np.array([1.0/s for s in grid_step] + [0.0], dtype=np.float32),
            'grid_p0':      np.array([*grid_p0   ,0.0],                  dtype=np.float32),
            'GFFParams':    np.array([r_damp, alpha_morse, 0.0, 0.0],    dtype=np.float32),
            'nstep':        np.int32(self.nstep),
        })
        
        # 3. Create buffers BEFORE generating kernel args
        if use_texture:
            print(f"MolecularDynamics::initGridFF() use_texture=True grid_shape={grid_shape} bspline_data.shape={bspline_data.shape} bspline_data.dtype={bspline_data.dtype}")
            print("Original bspline_data dimensions:", bspline_data.ndim)

            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = (40,40,200)   # better 
            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = ( grid_shape[2], grid_shape[1], grid_shape[0] )   # better 
            bspline_data = bspline_data.transpose(2,1,0,3).copy(); grid_shape = ( grid_shape[0], grid_shape[1], grid_shape[2] )   # better 
            print("!!!!! grid_shape: ", grid_shape)

            fmt = cl.ImageFormat(cl.channel_order.RGBA, cl.channel_type.FLOAT)
            tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, fmt, shape=tuple(grid_shape))
            
            #tex = cl.Image(self.ctx, cl.mem_flags.READ_ONLY, fmt, shape=tuple(grid_shape[::-1]))

            cl.enqueue_copy(self.queue, tex, bspline_data,  origin=(0, 0, 0), region=grid_shape)
            #cl.enqueue_copy(self.queue, tex, bspline_data,  origin=(0, 0, 0), region=bspline_data.shape[:3])
            self.buffer_dict['BsplinePLQH_tex'] = tex
            if bKernels:
                self.kernel_args_getNonBond_GridFF_Bspline_tex = self.generate_kernel_args("getNonBond_GridFF_Bspline_tex", bPrint=False)            
        else:
            print("MolecularDynamics::initGridFF() use_texture=False")

            #bspline_data = bspline_data.transpose(2,1,0,3).copy(); # grid_shape = ( grid_shape[2], grid_shape[1], grid_shape[0] )   # better 

            print("!!!!! grid_shape: ", grid_shape)

            buf = cl.Buffer(self.ctx, cl.mem_flags.READ_ONLY|cl.mem_flags.COPY_HOST_PTR, hostbuf=bspline_data)
            self.buffer_dict['BsplinePLQH'] = buf
            if bKernels:
                self.kernel_args_getNonBond_GridFF_Bspline = self.generate_kernel_args("getNonBond_GridFF_Bspline", bPrint=False)
        print("MolecularDynamics::initGridFF() DONE")

    def scan_1D( self, nsteps=100, d=[0.1,0.0,0.0], p0=[0.0,0.0,0.0], use_texture=False, mmff=None ):
        if mmff is None: mmff = self.mmff_list[0]
        # Initialize position and force grid_dataays
        pos    = np.zeros((nsteps,4), dtype=np.float32)
        forces = np.zeros((nsteps,4), dtype=np.float32)
        print("Running 1D force scan...")
        d  = np.array(d, dtype=np.float32)
        p0 = np.array(p0,  dtype=np.float32)
        for i in range(nsteps):
            # Update atom position along x-axis
            mmff.apos[0,:3] = p0 + d*i
            self.upload_all_systems()
            if use_texture:
                self.run_getNonBond_GridFF_Bspline_tex()
            else:
                self.run_getNonBond_GridFF_Bspline()
            pos_i, force_i = self.download_results()
            forces[i] = force_i[0]  # Get force on first (and only) atom
            pos[i] = pos_i[0]  # Store position for reference
            #if i % 10 == 0:  print(f"Step {i+1}/{nsteps}: x = {x[i]:.2f} Å, F = ({forces[i,0]:.3f}, {forces[i,1]:.3f}, {forces[i,2]:.3f}) kJ/mol/Å")
        return pos, forces

    def scan_2D( self, ns=(50,50), du=[0.1,0.0,0.0], dv=[0.0,0.1,0.0],  p0=[0.0,0.0,0.0], use_texture=False, mmff=None ):
        if mmff is None: mmff = self.mmff_list[0]
        pos    = np.zeros((ns[0],ns[1],4), dtype=np.float32)
        forces = np.zeros((ns[0],ns[1],4), dtype=np.float32)
        du = np.array(du, dtype=np.float32)
        dv = np.array(dv, dtype=np.float32)
        p0 = np.array(p0,  dtype=np.float32)
        print("Running 2D force scan...")
        for ix in range(ns[0]):
            for iy in range(ns[1]):
                p = p0 + du*ix + dv*iy
                mmff.apos[0,:3] = p
                self.upload_all_systems()
                if use_texture:
                    self.run_getNonBond_GridFF_Bspline_tex()
                else:
                    self.run_getNonBond_GridFF_Bspline()
                pos_i, force_i = self.download_results()
                forces[ix, iy,:] = force_i[0,: ]  # Get force on first (and only) atom
                pos   [ix, iy,:] = pos_i  [0,: ]  # Store position for reference
        return pos, forces

    def realloc_scan(self, n, na=-1):
        sz_f  = 4
        buffs = {
            "poss":    (sz_f*4 * n),
            "forces":  (sz_f*4 * n),
        }
        if na>0:
            buffs.update({
                "apos":   (sz_f*4 * na),
                "aREQs":  (sz_f*4 * na),
            })
        self.try_make_buffers(buffs)

    def scanNonBond(self, pos, force, REQH, ffpar, bRealloc=True ):
        n = len(pos)
        if bRealloc:  self.realloc_scan(n)
        # Upload data to GPU
        self.toGPU_( self.poss_buff,   pos)
        self.toGPU_( self.forces_buff, force)
        # Get kernel
        kernel = self.prg.scanNonBond
        # Set arguments
        kernel.set_args(
            np.int32(n),
            cl_array.vec.make_float4(*REQH),
            self.poss_buff,
            self.forces_buff,
            cl_array.vec.make_float8(*ffpar)
        )
        # Run kernel
        global_size = (n,)
        local_size = None
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        result = self.fromGPU_( self.forces_buff, shape=(n, 4))
        self.queue.finish()
        return result

    def scanNonBond2(self, pos, force, apos, aREQs, REQH0, ffpar, bRealloc=True, nPBC=None, lvec=None, name=""):
        n  = len(pos)
        na = len(apos)
        if bRealloc:  self.realloc_scan(n, na=na)
        # Upload data to GPU
        self.toGPU_( self.poss_buff,   pos)
        self.toGPU_( self.forces_buff, force)
        self.toGPU_( self.apos_buff,   apos)
        self.toGPU_( self.aREQs_buff,  aREQs)  
        # Get kernel

        if nPBC is not None:
            # Convert lvec to cl_Mat3 structure (3 float4 vectors)
            lvec_cl = np.zeros((3,4), dtype=np.float32)
            lvec_cl[:,:3] = lvec
            npbc_cl = np.zeros((4), dtype=np.int32)
            npbc_cl[:3] = nPBC
            #kernel = self.prg.scanNonBond2PBC
            kernel = self.prg.scanNonBond2PBC_2
            kernel.set_args(
                np.int32(n),
                cl_array.vec.make_float4(*REQH0),
                self.poss_buff,
                self.forces_buff,
                np.int32(na),
                self.apos_buff,
                self.aREQs_buff,
                cl_array.vec.make_float8(*ffpar),
                lvec_cl,
                npbc_cl,
            )
        else:
            kernel = self.prg.scanNonBond2        
            kernel.set_args(
                np.int32(n),
                cl_array.vec.make_float4(*REQH0),
                self.poss_buff,
                self.forces_buff,
                np.int32(na),
                self.apos_buff,
                self.aREQs_buff,
                cl_array.vec.make_float8(*ffpar)
            )        
        nloc=32
        local_size  = (nloc,)
        global_size = (clu.roundup_global_size(n, nloc),)
        T0 = time.time()
        cl.enqueue_nd_range_kernel(self.queue, kernel, global_size, local_size)
        self.queue.finish()
        T = time.time() - T0

        if nPBC is not None:
            npbc=(nPBC[0]*2+1)*(nPBC[1]*2+1)*(nPBC[2]*2+1)
            ntot = n*na*npbc
            #print(f"scanNonBond2() {name:<15} | {(T*1.e+9/ntot):>2.6f} [ns/op] {(T*1.e+12/ntot):>4.6f} [TOPS] | ntot: {ntot:<12}  np: {n:<6} na: {na:<6} nPBC({npbc:<6},{nPBC}) time: {T:3.6f} [s]")
            print(f"scanNonBond2PBC() {name:<15} | {T*1.e+9/ntot:>8.4f} [ns/op] {(ntot/(T*1.e+9)):>8.4f} [GOPS] | ntot: {ntot:>12} np: {n:>6} na: {na:>6} nPBC({npbc:>6},{nPBC}) time: {T:>8.4f} [s]")
        else:
            ntot = n*na
            #print(f"scanNonBond2() {name:<15} | {(T*1.e+9/ntot):>2.6f} [ns/op] {(T*1.e+12/ntot):>4.6f} [TOPS] | ntot: {ntot:<12}  np: {n:<6} na: {na:<6} nPBC({npbc:<6},{nPBC}) time: {T:3.6f} [s]")
            print(f"scanNonBond2() {name:<15} | {T*1.e+9/ntot:>8.4f} [ns/op] {(ntot/(T*1.e+9)):>8.4f} [GOPS] | ntot: {ntot:>12} np: {n:>6} na: {na:>6} time: {T:>8.4f} [s]")

        # Download results
        result = self.fromGPU_( self.forces_buff, shape=(n, 4))
        return result