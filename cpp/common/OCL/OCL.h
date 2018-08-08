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
    cl_mem_flags flags    = CL_MEM_READ_WRITE;
    char       * name      = NULL;
    // following is needed only for images
    bool img_dims = 0;
    int  nImg[3]  = {0,0,0};
    cl_image_format imageFormat;

    inline int initOnGPU ( cl_context& context ){
        int err;
        if( flags == CL_MEM_WRITE_ONLY ){  p_gpu = clCreateBuffer(context, flags, typesize * n, NULL,    &err);
        }else{                 if(p_cpu)  p_gpu = clCreateBuffer(context, flags, typesize * n, p_cpu,   &err); }
        return err;
    }

    inline int initOnGPUImage( cl_context& context ){
        int err;
        //p_gpu = clCreateBuffer(context, flags, typesize * n, NULL,  &err);
        p_gpu = clCreateImage2D(context, flags, &imageFormat, nImg[0],nImg[1], 0, p_cpu, &err);   // TODO: ??? nx=nImg[0] ny=nImg[1]  ???
        return err;
    }

    inline int setAsArg( cl_kernel& kernel, int i   ){ return clSetKernelArg(kernel, i, sizeof(cl_mem), &p_gpu );  };
    inline int fromGPU ( cl_command_queue& commands ){ return clEnqueueReadBuffer ( commands, p_gpu, CL_TRUE, 0, typesize * n, p_cpu, 0, NULL, NULL);  }
    inline int toGPU   ( cl_command_queue& commands ){ return clEnqueueWriteBuffer( commands, p_gpu, CL_TRUE, 0, typesize * n, p_cpu, 0, NULL, NULL ); }

    //inline setImageParams(  );

    inline OCLBuffer(){};
    inline OCLBuffer( char* name_, size_t n_, size_t typesize_, void * p_cpu_, cl_mem_flags flags_ ) :n(n_),typesize(typesize_),p_cpu(p_cpu_),flags(flags_),name(name_){};
};

class OCLsystem{
    // http://stackoverflow.com/questions/20105566/advantages-of-a-program-containing-several-opencl-kernels-versus-several-program
    public:
    cl_int           err;           // error code returned from OpenCL calls
    cl_device_id     device   = 0;        // compute device id
    cl_context       context  = 0;       // compute context
    cl_command_queue commands = 0;      // compute command queue
    cl_program       program  = 0;       // compute program - TODO FIXME: There could be more than one !!!!

    std::vector<cl_kernel> kernels;
    std::vector<OCLBuffer> buffers;

    void check_programSet (){ if(program ==0){ printf("ERROR OCLsystem program  not set \n"); exit(-1); } }
    void check_contextSet (){ if(context ==0){ printf("ERROR OCLsystem context  not set \n"); exit(-1); } }
    void check_deviceSet  (){ if(device  ==0){ printf("ERROR OCLsystem device   not set \n"); exit(-1); } }
    void check_commandsSet(){ if(commands==0){ printf("ERROR OCLsystem commands not set \n"); exit(-1); } }

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
        context  = clCreateContext(0, 1, &device, NULL, NULL, &err);  OCL_checkError(err, "Creating context");
        commands = clCreateCommandQueue(context, device, 0, &err);    OCL_checkError(err, "Creating command queue");
        return err;
    }

    int newKernel( char * name ){
        check_programSet();
        int err; kernels.push_back( clCreateKernel( program, name, &err ) );  OCL_checkError(err, "clCreateKernel"); return kernels.size()-1;
    }

    int newBuffer( char* name, size_t n, size_t typesize, void * p_cpu, cl_mem_flags flags ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, n, typesize, p_cpu, flags ) ); int i=buffers.size()-1; int err=buffers[i].initOnGPU(context); OCL_checkError(err, "initOnGPU"); return i;
    }

    int newBufferImage2D( char* name, size_t nx, size_t ny, size_t typesize, void * p_cpu, cl_mem_flags flags, cl_image_format imageFormat ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, nx*ny, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffers[i].img_dims    = 2;
        buffers[i].nImg[0]     = nx;
        buffers[i].nImg[1]     = ny;
        buffers[i].imageFormat = imageFormat;
        int err=buffers[i].initOnGPUImage(context); OCL_checkError(err, "initOnGPUImage");
        return i;
    }


    int initBuffers   (){ int err = CL_SUCCESS; for(int i=0; i<buffers.size(); i++){  err |= buffers[i].initOnGPU ( context );     }; return err; }
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

    // TODO : newProgram instead ?
    int buildProgram( char * fname ){
        check_deviceSet();
        char * kernelsource = getKernelSource( fname );
        // Create the comput program from the source buffer
        program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
        char tmpstr[1024];
        sprintf(tmpstr,"Creating program with %s", fname);
        OCL_checkError(err, tmpstr);
        free(kernelsource);
        err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
        if (err != CL_SUCCESS){
            printf( " ERROR in clBuildProgram %s \n", fname);
            OCL_buildProgramFailure( program, device );
            return -1;
        }
        //delete [] kernelsource; // TODO ??????
        return err;
    }

    inline int upload  (int i){ return buffers[i].toGPU(commands);    };
    inline int download(int i){ return buffers[i].fromGPU(commands);  };
    
    inline int upload  (int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].toGPU(commands);    };
    inline int download(int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].fromGPU(commands);  };
    
    inline int copy    (int from, int to, int from0, int to0, int n){ return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,from0,to0,n,0,NULL,NULL); };
    inline int copyBuff(int from, int to                           ){ int n=buffers[from].n; int n_=buffers[to].n; if(n_<n)n=n_; return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,0,0,n,0,NULL,NULL); };

    int download(){
        int err = CL_SUCCESS;
        for(int i=0; i<buffers.size(); i++ ){
            if( buffers[i].read_on_finish ){
                //printf("finish : reading buff %i \n", i);
                err |= buffers[i].fromGPU( commands );
            }
        }
        return err;
    }

    int finish(){
        int err;
        err = clFinish(commands);   OCL_checkError(err, "Waiting for kernel to finish");
        err |= download();
        return err;
    }

    void destroy(){
        clReleaseProgram(program);
        for(int i=0; i<kernels.size(); i++){ clReleaseKernel(kernels[i]);       }
        for(int i=0; i<buffers.size(); i++){ clReleaseMemObject(buffers[i].p_gpu); }
        //clReleaseKernel(kernel);
        clReleaseCommandQueue(commands);
        clReleaseContext(context);
    }
};

#define OCL_BUFF   1
#define OCL_INT    2
#define OCL_FLOAT  3
#define OCL_LBUFF  4

#define BUFFarg(X)   OCLarg( X, OCL_BUFF  )
#define INTarg(X)    OCLarg( X, OCL_INT   )
#define FLOATarg(X)  OCLarg( (float)X     )
#define LBUFFarg(X)  OCLarg( X, OCL_LBUFF )

class OCLarg{
    public:
    int  kind=0;
    union{
        float  f;
        int    i;
    };
    inline void setFloat(  float f_ ){ f=f_; kind=OCL_FLOAT; }
    inline void setInt  (  int   i_ ){ i=i_; kind=OCL_FLOAT; }
    inline void setBuff (  int   i_ ){ i=i_; kind=OCL_BUFF;  }
    OCLarg(){};
    OCLarg(  float f_            ):f(f_){ kind=OCL_FLOAT; }
    OCLarg(  int   i_, int kind_ ):i(i_), kind(kind_){ }
};

class OCLtask{
    public:
    OCLsystem  * cl;
    //cl_kernel  * kernel;
    size_t      ikernel   = 0;
    size_t      dim       = 1;
    size_t      global[3] = {0,0,0};
    size_t      local [3] = {0,0,0};

    std::vector<OCLarg> args;

    int useArgs(){
        int err = CL_SUCCESS;
        cl_kernel kernel = cl->kernels[ikernel];
        for(int i=0; i<args.size(); i++){
            OCLarg& arg = args[i];
            switch(arg.kind){
                case OCL_BUFF:  err |= clSetKernelArg( kernel, i, sizeof(cl_mem), &(cl->buffers[arg.i].p_gpu) );  OCL_checkError(err, "setAsArg"); break;
                //case OCL_BUFF:  err |= cl->buffers[arg.i].setAsArg( kernel, i );                                OCL_checkError(err, "setAsArg"); break;
                case OCL_INT:   err |= clSetKernelArg( kernel, i, sizeof(int),    &(arg.i) );                     OCL_checkError(err, "setAsArg"); break;
                case OCL_FLOAT: err |= clSetKernelArg( kernel, i, sizeof(float),  &(arg.f) );                     OCL_checkError(err, "setAsArg"); break;
                case OCL_LBUFF: err |= clSetKernelArg( kernel, i, arg.i,  NULL );                                 OCL_checkError(err, "setAsArg"); break;
            }
        }
        return err;
    }

    inline int enque_raw(  ){
        //printf("enque_raw %i %i (%i,%i) (%i,%i)\n", ikernel, dim, global[0],global[1], local[0],local[1]);
        if(local[0]==0){ return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, NULL,  0, NULL, NULL );   }
        else{            return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL );   }
    }

    virtual int enque( ){
        int err;
        if( args.size() > 0 ) useArgs();
        err = enque_raw( );  OCL_checkError(err, "enque_raw");
        return err;
    }

    void print_arg_list(){
        printf("kernel( ");
        for(int i=0; i<args.size(); i++){
            switch(args[i].kind){
                case OCL_INT:   printf( "int %i, ",      args[i].i ); break;
                case OCL_FLOAT: printf( "float %g, ",    args[i].f ); break;
                case OCL_BUFF:  printf( "buff[%i]:%s, ", args[i].i, cl->buffers[args[i].i].name );   break;
            }
        }
        printf(")\n");
    }

    inline void setup( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ){ cl=cl_; ikernel=ikernel_; dim=dim_;  global[0]=global_;global[1]=global_; local[0]=local_;local[1]=local_; };
    OCLtask          ( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ): cl(cl_),ikernel(ikernel_),dim(dim_){ global[0]=global_;global[1]=global_; local[0]=local_;local[1]=local_; };
    OCLtask(){};
};

// ========== Helper functions for coverting buffers to OpenCL format




#endif
