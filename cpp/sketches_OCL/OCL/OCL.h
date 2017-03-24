#ifndef  clUtils_h
#define  clUtils_h

#include <vector>
#include <CL/cl.h>

#include "OCLerrors.h"
#include "OCL_device_picker.h"

class OCLBuffer{
    public:
    void       * p_cpu = NULL;
    cl_mem       p_gpu;
    size_t       n        = 0;
    size_t       typesize = 0;
    bool         read_on_finish = false;
    cl_mem_flags flags = CL_MEM_READ_WRITE;

    inline cl_int initOnGPU ( cl_context& context ){
        cl_int err;
        if( flags == CL_MEM_WRITE_ONLY ){  p_gpu = clCreateBuffer(context, flags, typesize * n, NULL,    &err);
        }else{                  if(p_cpu)  p_gpu = clCreateBuffer(context, flags, typesize * n, p_cpu,   &err); }
        return err;
    }

    inline cl_int setAsArg( cl_kernel& kernel, int i  ){ return clSetKernelArg(kernel, i, sizeof(cl_mem), p_gpu );  };
    inline cl_int fromGPU ( cl_command_queue& commands ){ return clEnqueueReadBuffer( commands, p_gpu, CL_TRUE, 0, typesize * n, p_cpu, 0, NULL, NULL); }
};

class OCLArg{
    public:
    int kind      = 0;
    OCLBuffer * buff;
};

class OCLKernel{
    public:
    cl_kernel   kernel;
    cl_uint     dim       = 1;
    size_t      global[2] = {0,0};
    size_t      local [2] = {0,0};
    cl_uint     blocksize = 0;

    int nargs;
    OCLArg * args = NULL;

    //std::vector<OCLBuffer> buffers;
    /*
    inline cl_int setArgs (){
        cl_int err = 0;
        for(int i=0; i<nargs; i++){
            if( ){
            err | = clSetKernelArg(kernel, args, sizeof(cl_mem), &d_a);
            switch( args[i].kind ){
                case 0:
                case 1:

            err | = clSetKernelArg(kernel, args, sizeof(cl_mem), &d_a);
            }
        }
        return err;
    }
    */

    inline cl_int enque( cl_command_queue& commands ){
        return clEnqueueNDRangeKernel( commands, kernel, dim, NULL, global, local, 0, NULL, NULL );
    }
};

class OCLsystem{
    // http://stackoverflow.com/questions/20105566/advantages-of-a-program-containing-several-opencl-kernels-versus-several-program
    public:
    cl_int           err;           // error code returned from OpenCL calls
    cl_device_id     device;        // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    //cl_kernel        kernel;        // compute kernel

    std::vector<OCLKernel> kernels;
    std::vector<OCLBuffer> buffers;

    int init(){
        cl_uint deviceIndex = 0;
        //parseArguments(argc, argv, &deviceIndex);
        cl_device_id devices[MAX_DEVICES];
        unsigned numDevices = getDeviceList(devices);
        if (deviceIndex >= numDevices){  printf("Invalid device index (try '--list')\n"); return -1; }
        device = devices[deviceIndex];
        char name[MAX_INFO_STRING];
        getDeviceName(device, name);
        printf("\nUsing OpenCL device: %s\n", name);
        context = clCreateContext(0, 1, &device, NULL, NULL, &err);  OCL_checkError(err, "Creating context");
        commands = clCreateCommandQueue(context, device, 0, &err);   OCL_checkError(err, "Creating command queue");
        return err;
    }

    int initBuffers   (){ for(int i=0; i<buffers.size(); i++){ buffers[i].initOnGPU ( context );     } }
    //int releaseBuffers(){ for(int i=0; i<buffers; i++){ clReleaseMemObject(buffers[i].p_gpu); } }

    char * getKernelSource(char *filename){
        FILE *file = fopen(filename, "r");
        if (!file){ fprintf(stderr, "Error: Could not open kernel source file\n"); exit(-1); }
        fseek(file, 0, SEEK_END);
        int len = ftell(file) + 1;
        rewind(file);
        char *source = (char *)calloc(sizeof(char), len);
        if (!source){ fprintf(stderr, "Error: Could not allocate memory for source string\n"); exit(-1); }
        fread(source, sizeof(char), len, file);
        fclose(file);
        return source;
    }

    int buildProgram( char * fname ){
        char * kernelsource = getKernelSource( fname);
        // Create the comput program from the source buffer
        cl_program  program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
        char tmpstr[1024];
        sprintf(tmpstr,"Creating program with %s", fname);
        OCL_checkError(err, tmpstr);
        free(kernelsource);
        // Build the program
        err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
        if (err != CL_SUCCESS){
            size_t len;
            //char buffer[2048];
            printf("Error: Failed to build program executable!\n%s\n", OCL_err_code(err));
            clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(tmpstr), tmpstr, &len);
            printf("%s\n", tmpstr);
            return -1;
        }
        // Create the compute kernel from the program
        //kernel = clCreateKernel(program, name, &err);
        //OCL_checkError(err, tmpstr);
        return err;
    }

    /*
    int buildKernel( char * fname, int nnames, char * name ){
        cl_kernel kernel;
        buildKernel( kernel, fname, name );
        kernels[kernels.size()].kernel = kernel;
        //kernels.push_back(kernel);
    }
    */

    int finish(){
        err = clFinish(commands);
        OCL_checkError(err, "Waiting for kernel to finish");
        for(int i=0; i<buffers.size(); i++ ){
            if( buffers[i].read_on_finish ) buffers[i].fromGPU( commands );
        }
        return err;
    }

    int runKernelSequence( int n, int * ikerns ){
        for(int i=0; i<n; i++){
            err |= kernels[ikerns[i]].enque(commands);
            //OCL_checkError(err, "kernel ");
        }
        return err;
    }

    void destroy(){
        clReleaseProgram(program);
        for(int i=0; i<kernels.size(); i++){ clReleaseKernel(kernels[i].kernel);       }
        for(int i=0; i<buffers.size(); i++){ clReleaseMemObject(buffers[i].p_gpu); }
        //clReleaseKernel(kernel);
        clReleaseCommandQueue(commands);
        clReleaseContext(context);
    }
};

#endif
