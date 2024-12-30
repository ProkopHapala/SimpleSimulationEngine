#ifndef  OCL_Orb_cpp_h
#define  OCL_Orb_cpp_h

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_MINIMUM_OPENCL_VERSION 200

#include <CL/opencl.hpp>
#include "datatypes.h"
#include "datatypes2cl.h"
#include "Vec3.h"
#include "OrbSim_f.h"

#include <memory>
#include <string>
#include <stdexcept>

inline size_t roundUp(size_t size, size_t alignment) {
    return ((size + alignment - 1) / alignment) * alignment;
}

class OCL_Orb_cpp : public OrbSim_f {
public:
    // OpenCL objects
    cl::Context context;
    cl::CommandQueue queue;
    cl::Program program1;
    
    // Kernels
    cl::Kernel ker_dot_mat_vec_loc;
    cl::Kernel ker_dot_mat_vec_sparse;
    cl::Kernel ker_getTrussForces;
    cl::Kernel ker_updateJacobi_neighs;
    cl::Kernel ker_updateJacobi_mix;
    cl::Kernel ker_PD_perdictor;
    cl::Kernel ker_PD_corrector;
    
    // Buffers
    struct Buffers {
        // Points and properties
        cl::Buffer points;     // atoms state
        cl::Buffer ps1;        // working copy for solver
        cl::Buffer ps2;        // working copy for solver
        cl::Buffer dps;        // momentum in solver
        cl::Buffer forces;
        cl::Buffer vels;
        
        // Neighbor data
        cl::Buffer neighs;     // neighbor indices
        cl::Buffer neighBs;    // bond neighbors
        cl::Buffer neighB2s;   // secondary neighbors
        cl::Buffer params;     // neighbor parameters
        
        // Bond data
        cl::Buffer bparams;    // bond parameters
        cl::Buffer bforces;    // bond forces
        cl::Buffer bonds;      // bond connectivity
        
        // CG solver
        cl::Buffer Amat;       // matrix for CG
        cl::Buffer xvec;       // solution vector
        cl::Buffer yvec;       // temporary vector
    } buf;

    // Parameters
    int4   nDOFs    {0,0,0,0};       // number of DOFs (nPoints,nNeighMax,0,0)
    Quat4f MDpars   {1e-3,1e-4,0,0}; // MD parameters (dt, damping, cv, cf)
    Vec2f  bmix     {1.0,0.7};       // momentum-mixing parameters

    int initCLBuffs() {
        printf("initCLBuffs() nPoint %i nNeighMax %i \n", nPoint, nNeighMax);
        try {
            // Points and their properties
            buf.points  = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            buf.ps1     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            buf.ps2     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            buf.dps     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            buf.forces  = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            buf.vels    = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));
            
            // Neighbor-related buffers
            buf.neighs  = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int));
            buf.neighBs = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int2));
            buf.neighB2s= cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int));
            buf.params  = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(Quat4f));
            
            // Bond-related buffers
            buf.bparams = cl::Buffer(context, CL_MEM_READ_ONLY,  nBonds    * sizeof(Quat4f));
            buf.bforces = cl::Buffer(context, CL_MEM_READ_WRITE, nBonds    * sizeof(Quat4f));
            buf.bonds   = cl::Buffer(context, CL_MEM_READ_ONLY,  nBonds    * sizeof(int2));
            
            return 0;
        } catch (cl::Error& e) {
            printf("OpenCL error in initCLBuffs: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    int initCLBuffs_CG() {
        printf("initCLBuffs_CG() nPoint %i nNeighMax %i \n", nPoint, nNeighMax);
        try {
            buf.Amat = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint * nPoint * sizeof(float));
            buf.xvec = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint * sizeof(Quat4f));
            buf.yvec = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint * sizeof(Quat4f));
            return 0;
        } catch (cl::Error& e) {
            printf("OpenCL error in initCLBuffs_CG: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    void run_ocl(int niter, int upload_mask = 0b001, int download_mask = 0b001) {
        try {
            MDpars = Quat4f{dt, 1-damping, cv, cf};

            // Upload data to GPU
            if (upload_mask & 0b001) queue.enqueueWriteBuffer(buf.points, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
            if (upload_mask & 0b010) queue.enqueueWriteBuffer(buf.vels,   CL_TRUE, 0, nPoint * sizeof(Quat4f), vel);
            if (upload_mask & 0b100) queue.enqueueWriteBuffer(buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);

            // Execute kernels
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            
            for(int itr = 0; itr < niter; itr++) {
                queue.enqueueNDRangeKernel(ker_getTrussForces, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
                queue.enqueueNDRangeKernel(ker_PD_perdictor,   cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
            }

            // Download results from GPU
            if (download_mask & 0b001) queue.enqueueReadBuffer(buf.points, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
            if (download_mask & 0b010) queue.enqueueReadBuffer(buf.vels,   CL_TRUE, 0, nPoint * sizeof(Quat4f), vel);
            if (download_mask & 0b100) queue.enqueueReadBuffer(buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);

            queue.finish();
        } catch (cl::Error& e) {
            printf("OpenCL error in run_ocl: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

protected:
    cl::Program buildProgram(const std::string& filename) {
        FILE* file = fopen(filename.c_str(), "r");
        if (!file) {
            printf("Could not open kernel source file: %s\n", filename.c_str());
            throw std::runtime_error("Failed to open kernel file");
        }
        
        fseek(file, 0, SEEK_END);
        size_t size = ftell(file);
        rewind(file);
        
        std::string source(size, '\0');
        if (fread(&source[0], 1, size, file) != size) {
            printf("Failed to read kernel source file: %s\n", filename.c_str());
            fclose(file);
            throw std::runtime_error("Failed to read kernel file");
        }
        fclose(file);
        
        try {
            cl::Program program(context, source);
            program.build("-I. -cl-std=CL2.0");
            return program;
        } catch (cl::Error& e) {
            printf("OpenCL build error: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }
};

#endif