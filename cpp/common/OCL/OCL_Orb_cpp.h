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

    OCL_Orb_cpp() = default;

    void init(const char* cl_src_dir) {
        try {
            // Get platform and device
            std::vector<cl::Platform> platforms;
            cl::Platform::get(&platforms);
            if(platforms.empty()) {
                throw std::runtime_error("No OpenCL platforms found");
            }
            
            std::vector<cl::Device> devices;
            platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);
            if(devices.empty()) {
                throw std::runtime_error("No OpenCL GPU devices found");
            }
            
            // Create context and command queue
            context = cl::Context(devices[0]);
            queue   = cl::CommandQueue(context, devices[0]);
            
            // Initialize kernels
            char srcpath[1024];
            sprintf(srcpath, "%s/orb.cl", cl_src_dir);
            printf("Initializing OpenCL kernels from: %s\n", srcpath);
            
            initKernels();
            
            printf("OpenCL initialization complete\n");
        } catch (cl::Error& e) {
            printf("OpenCL error in init: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    void initKernels() {
        try {
            // Build the program
            program1 = buildProgram("cpp/common/OCL/cl/orb.cl");
            
            // Create kernels
            ker_getTrussForces      = cl::Kernel(program1, "getTrussForces");
            ker_updateJacobi_neighs = cl::Kernel(program1, "updateJacobi_neighs");
            ker_updateJacobi_mix    = cl::Kernel(program1, "updateJacobi_mix");
            ker_PD_perdictor        = cl::Kernel(program1, "PD_perdictor");
            ker_PD_corrector        = cl::Kernel(program1, "PD_corrector");
            ker_dot_mat_vec_loc     = cl::Kernel(program1, "dot_mat_vec_loc");
            ker_dot_mat_vec_sparse  = cl::Kernel(program1, "dot_mat_vec_sparse");
            
            printf("Successfully initialized OpenCL kernels\n");
        } catch (cl::Error& e) {
            printf("OpenCL error in initKernels: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }
    
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

    void run_getTrussForces(cl::Buffer* ps, cl::Buffer* fs ) {
        try {
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            ker_getTrussForces.setArg(0, nPoint);
            ker_getTrussForces.setArg(1, *ps );
            ker_getTrussForces.setArg(2, *fs );
            queue.enqueueNDRangeKernel(ker_getTrussForces, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
        } catch (cl::Error& e) {
            printf("OpenCL error in run_getTrussForces: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    void run_updateJacobi_neighs(cl::Buffer* ps_in, cl::Buffer* ps_out) {
        try {
            float inv_dt2 = 1.0f / (dt * dt);
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            
            ker_updateJacobi_neighs.setArg(0, nPoint);
            ker_updateJacobi_neighs.setArg(1, nNeighMax);
            ker_updateJacobi_neighs.setArg(2, *ps_in  );
            ker_updateJacobi_neighs.setArg(3, *ps_out );
            ker_updateJacobi_neighs.setArg(4, buf.dps );
            ker_updateJacobi_neighs.setArg(5, buf.neighs);
            ker_updateJacobi_neighs.setArg(6, buf.params);
            ker_updateJacobi_neighs.setArg(7, inv_dt2);
            
            queue.enqueueNDRangeKernel(ker_updateJacobi_neighs, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
        } catch (cl::Error& e) {
            printf("OpenCL error in run_updateJacobi_neighs: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    void run_updateJacobi_mix( cl::Buffer* ps_in, cl::Buffer* ps_out) {
        try {
            float inv_dt2 = 1.0f / (dt * dt);
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            
            ker_updateJacobi_mix.setArg(0, nPoint);
            ker_updateJacobi_mix.setArg(1, nNeighMax);
            ker_updateJacobi_mix.setArg(2, *ps_in );
            ker_updateJacobi_mix.setArg(3, *ps_out );
            ker_updateJacobi_mix.setArg(4, buf.dps);
            ker_updateJacobi_mix.setArg(5, buf.neighs);
            ker_updateJacobi_mix.setArg(6, buf.params);
            ker_updateJacobi_mix.setArg(7, inv_dt2);
            ker_updateJacobi_mix.setArg(8, bmix);
            
            queue.enqueueNDRangeKernel(ker_updateJacobi_mix, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
        } catch (cl::Error& e) {
            printf("OpenCL error in run_updateJacobi_mix: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    double run_updateJacobi_smart( cl::Buffer* psa, cl::Buffer* psb, int itr) {
        const float bmix_end   = 0.75f;
        const float bmix_start = 0.55f;
        const int   nitr_start = 3;
        const int   nitr_end   = 10;
        
        if(itr < nitr_start) {
            run_updateJacobi_neighs(psa, psb);
        } else {
            if(itr > nitr_end) {
                bmix.y = bmix_end;
            } else {
                bmix.y = bmix_start + (bmix_end - bmix_start) * (itr - nitr_start) / (nitr_end - nitr_start);
            }
            run_updateJacobi_mix(psa, psb);
        }
        return bmix.y;
    }

    void run_PD_perdictor( cl::Buffer* ps_out ) {
        try {
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            
            ker_PD_perdictor.setArg(0, nPoint);
            ker_PD_perdictor.setArg(1, buf.points);
            ker_PD_perdictor.setArg(2, *ps_out );
            ker_PD_perdictor.setArg(3, buf.forces);
            ker_PD_perdictor.setArg(4, buf.vels);
            ker_PD_perdictor.setArg(5, dt);
            
            queue.enqueueNDRangeKernel(ker_PD_perdictor, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
        } catch (cl::Error& e) {
            printf("OpenCL error in run_PD_perdictor: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    void run_PD_corrector( cl::Buffer* ps_in, cl::Buffer* ps_out  ) {
        try {
            size_t local_work_size = 32;
            size_t global_work_size = roundUp(nPoint, local_work_size);
            ker_PD_corrector.setArg(0, nPoint     );
            ker_PD_corrector.setArg(1, *ps_in     );
            ker_PD_corrector.setArg(2, buf.points );
            ker_PD_corrector.setArg(3, buf.vels   );
            ker_PD_corrector.setArg(4, dt);
            queue.enqueueNDRangeKernel(ker_PD_corrector, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));
        } catch (cl::Error& e) {
            printf("OpenCL error in run_PD_corrector: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }



    void run_projective_dynamics(int nSolverIters, int ialg) {
        cl::Buffer* psa =  ps1;      
        cl::Buffer* psb =  ps2;       
        run_PD_perdictor(psa);
        switch(ialg) {
            case 0: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_neighs( psa, psb );   std::swap(psa, psb); } break;
            case 1: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_mix   ( psa, psb );   std::swap(psa, psb); } break;
            case 2: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_smart ( psa, psb, i); std::swap(psa, psb); } break;
        }
        run_PD_corrector(psa);
    }

    void run_PDcl(int niter, int nSolverIters, int upload_mask = 0b001, int download_mask = 0b001, int ialg = 2) {
        try {
            MDpars = Quat4f{dt, 1-damping, cv, cf};
            
            if(upload_mask & 0b001) queue.enqueueWriteBuffer(buf.points, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
            if(upload_mask & 0b010) queue.enqueueWriteBuffer(buf.vels,   CL_TRUE, 0, nPoint * sizeof(Quat4f), vel);
            if(upload_mask & 0b100) queue.enqueueWriteBuffer(buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);
            
            for(int itr = 0; itr < niter; itr++) {
                run_projective_dynamics(nSolverIters, ialg);
            }
            
            if(download_mask & 0b001) queue.enqueueReadBuffer(buf.points, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
            if(download_mask & 0b010) queue.enqueueReadBuffer(buf.vels,   CL_TRUE, 0, nPoint * sizeof(Quat4f), vel);
            if(download_mask & 0b100) queue.enqueueReadBuffer(buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);
            
            queue.finish();
        } catch (cl::Error& e) {
            printf("OpenCL error in run_PDcl: %s (%d)\n", e.what(), e.err());
            throw;
        }
    }

    Vec2f run_SolverConvergence(int nSolverIters, int ialg, bool bPrint = false) {
        try {
            int psa = 1, psb = 2;  // Using 1 for ps1 and 2 for ps2 buffers
            
            // Initialize forces to zero
            std::fill(forces, forces + nPoint, Quat4fZero);
            
            // Upload initial data
            queue.enqueueWriteBuffer(buf.ps1, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
            queue.enqueueWriteBuffer(buf.dps, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);
            
            Vec2f fe = Vec2fNAN;
            Vec2f fe1 = Vec2fNAN;
            
            for(int i = 0; i < nSolverIters; i++) {
                // Run solver step
                switch(ialg) {
                    case 0: run_updateJacobi_neighs(psa, psb); break;
                    case 1: run_updateJacobi_mix(psa, psb);    break;
                    case 2: run_updateJacobi_smart(psa, psb, i); break;
                }
                
                // Get forces
                run_getTrussForces(psb, 1);  // 1 for main forces buffer
                
                // Download results
                queue.enqueueReadBuffer((psb == 1) ? buf.ps1 : buf.ps2, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);
                queue.enqueueReadBuffer(buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);
                queue.finish();
                
                // Calculate energy and force
                Vec2f fe_ = Vec2fZero;
                for(int ia = 0; ia < nPoint; ia++) {
                    fe_.x += forces[ia].f.norm2();
                    fe_.y += forces[ia].e;
                }
                
                std::swap(psa, psb);
                if(i == 0) fe1 = fe_;
                fe = fe_;
            }
            
            if(bPrint) {
                printf("RESULT OCL_Orb::run_SolverConvergence()  ialg: %i bmix(%.2f,%.2f) nstep: %3i Estart: %.2e Eend: %.2e Eend/Estart: %g\n", 
                    ialg, bmix.x, bmix.y, nSolverIters, fe1.y, fe.y, fe.y/fe1.y);
            }
            
            return fe;
        } catch (cl::Error& e) {
            printf("OpenCL error in run_SolverConvergence: %s (%d)\n", e.what(), e.err());
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