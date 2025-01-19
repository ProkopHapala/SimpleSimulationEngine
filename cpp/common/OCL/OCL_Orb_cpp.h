#ifndef  OCL_Orb_cpp_h
#define  OCL_Orb_cpp_h

#define CL_HPP_TARGET_OPENCL_VERSION 200
#define CL_HPP_MINIMUM_OPENCL_VERSION 200

#ifdef DEBUG
    #define CL_CHECK(err) if(err != CL_SUCCESS) { printf("OpenCL error: %d\n", err); }
#else
    #define CL_CHECK(err)
#endif

#include <CL/opencl.hpp>
#include "datatypes.h"
#include "datatypes2cl.h"
#include "Vec3.h"
#include "TrussDynamics_f.h"

#include <memory>
#include <string>
#include <stdexcept>

inline size_t roundUp(size_t size, size_t alignment) {
    return ((size + alignment - 1) / alignment) * alignment;
}

class OCL_Orb_cpp : public TrussDynamics_f {
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
        std::vector<cl::Platform> platforms;
        cl_int err = cl::Platform::get(&platforms);                                CL_CHECK(err);
        if(platforms.empty()) {
            printf("No OpenCL platforms found\n");
            return;
        }
        
        std::vector<cl::Device> devices;
        err = platforms[0].getDevices(CL_DEVICE_TYPE_GPU, &devices);              CL_CHECK(err);
        if(devices.empty()) {
            printf("No OpenCL GPU devices found\n");
            return;
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
    }

    void initKernels() {
        // Build the program
        program1 = buildProgram("cpp/common/OCL/cl/orb.cl");
        
        // Create kernels
        cl_int err;
        ker_getTrussForces      = cl::Kernel(program1, "getTrussForces");         CL_CHECK(err);
        ker_updateJacobi_neighs = cl::Kernel(program1, "updateJacobi_neighs");    CL_CHECK(err);
        ker_updateJacobi_mix    = cl::Kernel(program1, "updateJacobi_mix");       CL_CHECK(err);
        ker_PD_perdictor        = cl::Kernel(program1, "PD_perdictor");           CL_CHECK(err);
        ker_PD_corrector        = cl::Kernel(program1, "PD_corrector");           CL_CHECK(err);
        ker_dot_mat_vec_loc     = cl::Kernel(program1, "dot_mat_vec_loc");        CL_CHECK(err);
        ker_dot_mat_vec_sparse  = cl::Kernel(program1, "dot_mat_vec_sparse");     CL_CHECK(err);
        
        printf("Successfully initialized OpenCL kernels\n");
    }
    
    int initCLBuffs() {
        printf("initCLBuffs() nPoint %i nNeighMax %i \n", nPoint, nNeighMax);
        cl_int err;
        
        // Points and their properties
        buf.points  = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        buf.ps1     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        buf.ps2     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        buf.dps     = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        buf.forces  = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        buf.vels    = cl::Buffer(context, CL_MEM_READ_WRITE, nPoint    * sizeof(Quat4f));  CL_CHECK(err);
        
        // Neighbor-related buffers
        buf.neighs  = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int));     CL_CHECK(err);
        buf.neighBs = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int2));    CL_CHECK(err);
        buf.neighB2s= cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(int));     CL_CHECK(err);
        buf.params  = cl::Buffer(context, CL_MEM_READ_ONLY,  nNeighTot * sizeof(Quat4f));  CL_CHECK(err);
        
        // Bond-related buffers
        buf.bparams = cl::Buffer(context, CL_MEM_READ_ONLY,  nBonds    * sizeof(Quat4f));  CL_CHECK(err);
        buf.bforces = cl::Buffer(context, CL_MEM_READ_WRITE, nBonds    * sizeof(Quat4f));  CL_CHECK(err);
        buf.bonds   = cl::Buffer(context, CL_MEM_READ_ONLY,  nBonds    * sizeof(int2));    CL_CHECK(err);
        
        return 0;
    }

    void run_getTrussForces(cl::Buffer* ps, cl::Buffer* fs) {
        size_t local_work_size = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        cl_int err=0;
        cl::Kernel& ker = ker_getTrussForces;
        err = ker.setArg(0, nPoint);                               CL_CHECK(err);
        err = ker.setArg(1, *ps);                                  CL_CHECK(err);
        err = ker.setArg(2, *fs);                                  CL_CHECK(err);
        err = queue.enqueueNDRangeKernel(ker, cl::NullRange,  cl::NDRange(global_work_size), cl::NDRange(local_work_size));         CL_CHECK(err);
    }

    void run_updateJacobi_neighs(cl::Buffer* ps_in, cl::Buffer* ps_out) {
        float inv_dt2 = 1.0f / (dt * dt);
        size_t local_work_size = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        cl_int err=0;
        cl::Kernel& ker = ker_updateJacobi_neighs;
        
        err = ker.setArg(0, nPoint);                         CL_CHECK(err);
        err = ker.setArg(1, nNeighMax);                      CL_CHECK(err);
        err = ker.setArg(2, *ps_in);                         CL_CHECK(err);
        err = ker.setArg(3, *ps_out);                        CL_CHECK(err);
        err = ker.setArg(4, buf.dps);                        CL_CHECK(err);
        err = ker.setArg(5, buf.neighs);                     CL_CHECK(err);
        err = ker.setArg(6, buf.params);                     CL_CHECK(err);
        err = ker.setArg(7, inv_dt2);                        CL_CHECK(err);
        err = queue.enqueueNDRangeKernel(ker, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));        CL_CHECK(err);
    }

    void run_updateJacobi_mix(cl::Buffer* ps_in, cl::Buffer* ps_out) {
        float inv_dt2 = 1.0f / (dt * dt);
        size_t local_work_size = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        cl::Kernel& ker = ker_updateJacobi_mix;
        cl_int err;
        err = ker.setArg(0, nPoint);                            CL_CHECK(err);
        err = ker.setArg(1, nNeighMax);                         CL_CHECK(err);
        err = ker.setArg(2, *ps_in);                            CL_CHECK(err);
        err = ker.setArg(3, *ps_out);                           CL_CHECK(err);
        err = ker.setArg(4, buf.dps);                           CL_CHECK(err);
        err = ker.setArg(5, buf.neighs);                        CL_CHECK(err);
        err = ker.setArg(6, buf.params);                        CL_CHECK(err);
        err = ker.setArg(7, inv_dt2);                           CL_CHECK(err);
        err = ker.setArg(8, bmix);                              CL_CHECK(err);
        err = queue.enqueueNDRangeKernel(ker, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));        CL_CHECK(err);
    }

    double run_updateJacobi_smart(cl::Buffer* psa, cl::Buffer* psb, int itr) {
        const float bmix_end   = 0.75f;
        const float bmix_start = 0.55f;
        const int   nitr_start = 3;
        const int   nitr_end   = 10;
        if(itr < nitr_start) {
            run_updateJacobi_neighs(psa, psb);
        } else {
            if(itr > nitr_end){ bmix.y = bmix_end;  } 
            else              { bmix.y = bmix_start + (bmix_end - bmix_start) * (itr - nitr_start) / (nitr_end - nitr_start); }
            run_updateJacobi_mix(psa, psb);
        }
        return bmix.y;
    }

    void run_PD_perdictor(cl::Buffer* ps_out) {
        size_t local_work_size = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        cl_int err;
        cl::Kernel& ker = ker_PD_perdictor;
        err = ker.setArg(0, nPoint);                                CL_CHECK(err);
        err = ker.setArg(1, buf.points);                            CL_CHECK(err);
        err = ker.setArg(2, *ps_out);                               CL_CHECK(err);
        err = ker.setArg(3, buf.forces);                            CL_CHECK(err);
        err = ker.setArg(4, buf.vels);                              CL_CHECK(err);
        err = ker.setArg(5, dt);                                    CL_CHECK(err);
        err = queue.enqueueNDRangeKernel(ker, cl::NullRange,  cl::NDRange(global_work_size), cl::NDRange(local_work_size));        CL_CHECK(err);
    }

    void run_PD_corrector(cl::Buffer* ps_in, cl::Buffer* ps_out) {
        size_t local_work_size = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        cl_int err;
        cl::Kernel& ker = ker_PD_corrector;
        err = ker.setArg(0, nPoint);                                CL_CHECK(err);
        err = ker.setArg(1, *ps_in);                                CL_CHECK(err);
        err = ker.setArg(2, buf.points);                            CL_CHECK(err);
        err = ker.setArg(3, buf.vels);                              CL_CHECK(err);
        err = ker.setArg(4, dt);                                    CL_CHECK(err);
        err = queue.enqueueNDRangeKernel(ker, cl::NullRange, cl::NDRange(global_work_size), cl::NDRange(local_work_size));        CL_CHECK(err);
    }

    void run_projective_dynamics(int nSolverIters, int ialg) {
        cl::Buffer* psa = &buf.ps1;      
        cl::Buffer* psb = &buf.ps2;       
        run_PD_perdictor(psa);
        switch(ialg) {
            case 0: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_neighs( psa, psb );   std::swap(psa, psb); } break;
            case 1: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_mix   ( psa, psb );   std::swap(psa, psb); } break;
            case 2: for(int i = 0; i < nSolverIters; i++) { run_updateJacobi_smart ( psa, psb, i); std::swap(psa, psb); } break;
        }
        run_PD_corrector(psa);
    }

    void run_PDcl(int niter, int nSolverIters, int upload_mask = 0b001, int download_mask = 0b001, int ialg = 2) {
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
    }

    Vec2f run_SolverConvergence(int nSolverIters, int ialg, bool bPrint = false) {
        cl::Buffer* psa =  &buf.ps1;      
        cl::Buffer* psb =  &buf.ps2;   
        // Initialize forces to zero
        std::fill(forces, forces + nPoint, Quat4fZero);
        cl_int err;
        err = queue.enqueueWriteBuffer(buf.ps1, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);  CL_CHECK(err);
        err = queue.enqueueWriteBuffer(buf.dps, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);  CL_CHECK(err);
        Vec2f fe  = Vec2fNAN;
        Vec2f fe1 = Vec2fNAN;
        for(int i = 0; i < nSolverIters; i++) {
            // Run solver step
            switch(ialg) {
                case 0: run_updateJacobi_neighs(psa, psb); break;
                case 1: run_updateJacobi_mix   (psa, psb);    break;
                case 2: run_updateJacobi_smart (psa, psb, i); break;
            }
            run_getTrussForces(psb, &buf.forces);  // 1 for main forces buffer
            err = queue.enqueueReadBuffer( *psb, CL_TRUE, 0, nPoint * sizeof(Quat4f), points);  CL_CHECK(err);
            err = queue.enqueueReadBuffer( buf.forces, CL_TRUE, 0, nPoint * sizeof(Quat4f), forces);  CL_CHECK(err);
            queue.finish();
            Vec2f fe_ = Vec2fZero;
            for(int ia = 0; ia < nPoint; ia++) { fe_.x += forces[ia].f.norm2(); fe_.y += forces[ia].e;  }
            std::swap(psa, psb);
            if(i == 0) fe1 = fe_;
            fe = fe_;
        }
        if(bPrint) {  printf("RESULT OCL_Orb::run_SolverConvergence()  ialg: %i bmix(%.2f,%.2f) nstep: %3i Estart: %.2e Eend: %.2e Eend/Estart: %g\n",  ialg, bmix.x, bmix.y, nSolverIters, fe1.y, fe.y, fe.y/fe1.y); }
        return fe;
    }

protected:
    cl::Program buildProgram(const std::string& filename) {
        FILE* file = fopen(filename.c_str(), "r");
        if (!file) {
            printf("Failed to open %s\n", filename.c_str());
            return cl::Program();
        }

        fseek(file, 0, SEEK_END);
        size_t length = ftell(file);
        fseek(file, 0, SEEK_SET);

        std::string source(length + 1, '\0');
        fread(&source[0], 1, length, file);
        fclose(file);

        cl::Program::Sources sources;
        sources.push_back({source.c_str(), source.length()});

        cl_int err;
        cl::Program program = cl::Program(context, sources, &err);                CL_CHECK(err);
        err = program.build({});                                                  CL_CHECK(err);
        
        if (err != CL_SUCCESS) {
            printf("Error building program: %s\n", program.getBuildInfo<CL_PROGRAM_BUILD_LOG>().c_str());
            return cl::Program();
        }

        return program;
    }
};

#endif