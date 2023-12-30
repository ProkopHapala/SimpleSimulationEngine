#ifndef  clUtils_h
#define  clUtils_h

#include <string>
#include <vector>
#include <unordered_map>

//#define CL_TARGET_OPENCL_VERSION 200

#include <CL/cl.h>

#include "OCLerrors.h"
#include "OCL_device_picker.h"

#include "datatypes.h"


//=======================================================================
//=======================================================================

class OCLBuffer{public:

    void       * p_cpu = NULL;
    cl_mem       p_gpu    = 0;
    size_t       n        = 0;
    size_t       typesize = 0;
    bool         read_on_finish = false;
    cl_mem_flags flags    = CL_MEM_READ_WRITE;
    //char       * name     = NULL;
    std::string   name;
    // following is needed only for images
    int  img_dims = 0;
    int  nImg[3]  = {0,0,0};
    cl_image_format imageFormat;

    inline size_t byteSize(){ return typesize * n; }

    inline int initOnGPU ( cl_context& context ){
        int err=0;
        //printf( "initOnGPU() buff_size %li | n %li typesize %li \n", byteSize(), n, typesize );
        if( (flags&CL_MEM_COPY_HOST_PTR)||(flags&CL_MEM_USE_HOST_PTR) ){
                p_gpu = clCreateBuffer(context, flags, byteSize(), p_cpu,   &err);
        }else{  p_gpu = clCreateBuffer(context, flags, byteSize(), 0,       &err); }
        //printf( "initOnGPU p_gpu: %li \n", (long)p_gpu );
        return err;
    }

    inline int initOnGPUImage( cl_context& context ){
        int err=0;
        //p_gpu = clCreateBuffer(context, flags, typesize * n, NULL,  &err);
        cl_image_desc img;
        img.image_array_size  = 0;
        img.image_row_pitch   = 0;
        img.image_slice_pitch = 0;
        img.num_mip_levels    = 0;
        img.num_samples       = 0;
        img.buffer            = 0; // cl_mem :  From Buffer ?

        switch(img_dims){
            case 2:
                printf( " initOnGPUImage: clCreateImage2D \n" );
                img.image_type  = CL_MEM_OBJECT_IMAGE2D;
                img.image_width  = nImg[0];
                img.image_height = nImg[1];
                img.image_depth  = 0;
                p_gpu = clCreateImage  (context, flags, &imageFormat, &img,  p_cpu, &err);
                //p_gpu = clCreateImage2D(context, flags, &imageFormat, nImg[0],nImg[1],          0,    p_cpu, &err);   // TODO: ??? nx=nImg[0] ny=nImg[1]  ???
                break;
            case 3:
                //printf( " initOnGPUImage: clCreateImage3D \n" );
                img.image_type  = CL_MEM_OBJECT_IMAGE3D;
                img.image_width  = nImg[0];
                img.image_height = nImg[1];
                img.image_depth  = nImg[2];
                p_gpu = clCreateImage  (context, flags, &imageFormat, &img,   p_cpu, &err);
                //p_gpu = clCreateImage3D(context, flags, &imageFormat, nImg[0],nImg[1],nImg[2], 0, 0, p_cpu, &err);   // TODO: ??? nx=nImg[0] ny=nImg[1]  ???
                //printf( "initOnGPUImage( flags %li, imageFormat{%i,%i} nImg(%i,%i,%i) \n", flags, imageFormat.image_channel_data_type, imageFormat.image_channel_order, nImg[0],nImg[1],nImg[2] );
                break;
        }
        //p_gpu = clCreateImage  (context, flags, &imageFormat, &img,                          p_cpu, &err);
        //printf( "initOnGPUImage img_dims: %i p_gpu: %li \n", img_dims, (long)p_gpu );
        return err;
    }

    inline int setAsArg( cl_kernel& kernel, int i   ){ return clSetKernelArg(kernel, i, sizeof(cl_mem), &p_gpu );  };

    inline int fromGPU ( cl_command_queue& commands,       void* cpu_data, int n_=-1, int i0=0 ){ 
        if(img_dims==0){
            if(n_<0)n_=n; 
            //printf( "fromGPU() nbutes=%i n_ %i i0 %i typesize %li p_gpu=%li cpu_data=%li \n", n_*typesize,  n_, i0, typesize, (long)p_gpu, (long)cpu_data ); // DEBUG
            //for(int i=0;i<n_*typesize;i++){ ((char*)cpu_data)[i]=0; } // DEBUG
            int iret = clEnqueueReadBuffer( commands, p_gpu, CL_TRUE, i0*typesize, n_*typesize, cpu_data, 0, NULL, NULL ); 
            //OCL_checkError(iret,"fromGPU()");
            return iret;
        }else{
            size_t offset[4]{0,0,0,0};
            size_t region[4]{nImg[0],nImg[1],nImg[2],0};
            return clEnqueueReadImage ( commands, p_gpu, CL_TRUE, offset, region, 0,0,  cpu_data, 0, NULL, NULL ); 
        }
    }

    inline int toGPU   ( cl_command_queue& commands, const void* cpu_data, int n_=-1, int i0=0 ){ 
        if(img_dims==0){
            if(n_<0)n_=n; 
            return clEnqueueWriteBuffer( commands, p_gpu, CL_TRUE, i0*typesize, typesize*n_, cpu_data, 0, NULL, NULL );
        }else{
            size_t offset[4]{0,0,0,0};
            size_t region[4]{nImg[0],nImg[1],nImg[2],0};
            return clEnqueueWriteImage ( commands, p_gpu, CL_TRUE, offset, region, 0,0,  cpu_data, 0, NULL, NULL ); 
        }
    }
    
    inline int fromGPU ( cl_command_queue& commands ){ return fromGPU ( commands, p_cpu ); }
    inline int toGPU   ( cl_command_queue& commands ){ return toGPU   ( commands, p_cpu ); }
    //inline setImageParams(  );

    inline OCLBuffer(){};
    inline OCLBuffer( const char* name_, size_t n_, size_t typesize_, void * p_cpu_, cl_mem_flags flags_=CL_MEM_READ_WRITE ) :n(n_),typesize(typesize_),p_cpu(p_cpu_),flags(flags_),name(name_){};
    inline int release(){ return clReleaseMemObject( p_gpu ); p_gpu=0; }
    //inline ~OCLBuffer(){ if(p_gpu) release(); };   // This makes some problem with de-allocation
}; // OCLBuffer

//=======================================================================
//=======================================================================

#define OCL_BUFF   1
#define OCL_INT    2
#define OCL_FLOAT  3
#define OCL_LBUFF  4
#define OCL_PTR    5

#define BUFFarg(X)  OCLarg( (int)X, OCL_BUFF  )
#define INTarg(X)   OCLarg( (int)X, OCL_INT   )
#define FLOATarg(X) OCLarg( (float)X     )
#define LBUFFarg(X) OCLarg( (int)X, OCL_LBUFF )
#define PTRarg(X,n) OCLarg( (void*)X, OCL_PTR, n )
#define REFarg(X)   OCLarg( (void*)&X, sizeof(X) )

#define _useArg(X)  OCLsystem::useArg_( (void*)&X, sizeof(X) )

class OCLarg{  public:
    int  kind=0;
    int  nbytes=0;
    union{
        float  f;
        int    i;
        void*  ptr;
    };
    inline void setFloat(  float f_ ){ f=f_; kind=OCL_FLOAT; }
    inline void setInt  (  int   i_ ){ i=i_; kind=OCL_INT  ; }
    inline void setBuff (  int   i_ ){ i=i_; kind=OCL_BUFF;  }
    OCLarg(){};
    OCLarg( float f_                ):f(f_), kind(OCL_FLOAT) {}//  printf( "!!!!! OCLarg(f %f kind %i )\n", f, kind ); }
    OCLarg( int   i_,   int kind_   ):i(i_), kind(kind_)     {}//  printf( "!!!!! OCLarg(i %i kind %i )\n", i, kind ); }
    OCLarg( void* ptr_, int nbytes_ ):ptr(ptr_), kind(OCL_PTR), nbytes(nbytes_){} //{  printf( "!!!!! OCLarg(ptr %li kind %i  nbytes %i )\n", (long)ptr, kind, nbytes ); }
}; // class OCLarg

//=======================================================================
//=======================================================================

class OCLsystem;

class OCLtask{ public:
    OCLsystem*  cl=0;
    //cl_kernel   kernel;
    size_t  ikernel  = 0;
    size_t  dim      = 1;
    //size_t  global[3] = {0,0,0};
    //size_t  local [3] = {0,0,0};
    size_t4  global;
    size_t4  local;
    std::vector<OCLarg> args;

    void roundSizes(){  
        size_t* global_=(size_t*)&global;
        size_t* local_ =(size_t*)&local;
        for(int i=0;i<dim;i++){
            int nL=local_[i]; if(nL<1){nL==1;}
            global_[i] = (((int)(global_[i]/nL))+1)*nL;
            printf("OCLtask roundSizes[%li/%i] \n", global_[i], nL );
        }
    }

    int  useArgs();
    void print_arg_list();
    int enque_raw();

    virtual int enque( ){
        int err=0;
        if( args.size() > 0 ) useArgs();
        err = enque_raw( );  OCL_checkError(err, "enque_raw");
        return err;
    }

    //inline void setup4( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t4 global_, size_t4 local_ ){ cl=cl_; ikernel=ikernel_; dim=dim_;  global[0]=global_.x;global[1]=global_.y;global[2]=global_.z; local[0]=local_.x;local[1]=local_.y;local[2]=local_.z; };
    //inline void setup( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ){ cl=cl_; ikernel=ikernel_; dim=dim_;  global[0]=global_;global[1]=global_;global[2]=global_; local[0]=local_;local[1]=local_;local[2]=local_; };
    //OCLtask          ( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t global_, size_t local_ ): cl(cl_),ikernel(ikernel_),dim(dim_){ global[0]=global_;global[1]=global_;global[2]=global_; local[0]=local_;local[1]=local_;local[2]=local_; };
    OCLtask()=default;
    OCLtask          ( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t4 global_, size_t4 local_ ):cl(cl_),ikernel(ikernel_),dim(dim_),global(global_),local(local_){};
    
    //OCLtask( &cl, iker, 2, n, 16    ); 
    OCLtask          ( OCLsystem  * cl_, size_t ikernel_, size_t dim_, int nglob, int nloc ):cl(cl_),ikernel(ikernel_),dim(dim_),global( (size_t4){(size_t)nglob,(size_t)nglob,(size_t)nglob,(size_t)nglob}),local((size_t4){(size_t)nloc,(size_t)nloc,(size_t)nloc,(size_t)nloc}){};
    
    inline void setup( OCLsystem  * cl_, size_t ikernel_, size_t dim_, size_t4 global_, size_t4 local_ ){cl=cl_; ikernel=ikernel_; dim=dim_; global=global_; local=local_; };

    //OCLtask(){}
}; // class OCLtask

//=======================================================================
//=======================================================================

class OCLsystem{ public:
    // http://stackoverflow.com/questions/20105566/advantages-of-a-program-containing-several-opencl-kernels-versus-several-program
    //cl_int           err;              // error code returned from OpenCL calls
    cl_device_id     device   = 0;      // compute device id
    cl_context       context  = 0;      // compute context
    cl_command_queue commands = 0;      // compute command queue
    cl_program       program  = 0;      // compute program - TODO FIXME: There could be more than one !!!!

    std::vector<cl_kernel> kernels;
    std::vector<OCLBuffer> buffers;
    std::vector<OCLtask*>  tasks;

    std::unordered_map<std::string,int> kernel_dict;
    std::unordered_map<std::string,int> buffer_dict;
    std::unordered_map<std::string,int> task_dict;

    cl_kernel current_kernel;
    OCLtask* currentTask=0;
    int argCounter;

    OCLtask* getTask(const char* name){ 
        auto it = task_dict.find(name);
        if(it != task_dict.end() ){ 
            return tasks[ it->second ];
        }
        printf( "ERROR in OCLsystem::getTask() task '%s' not found. Known tasks are: \n", name );
        for( auto it=task_dict.begin(); it!=task_dict.end(); it++ ){ printf( " %s \n", it->first.c_str() ); }
        exit(0);
        return 0;
    }    

    OCLtask* getOrMakeTask(const char* name){ 
        auto it = task_dict.find(name);
        if(it != task_dict.end() ){ 
            return tasks[ it->second ];
        }
        
        //printf( "OCL task %s not found \n", name );
        return 0;
    }

    int roundUpSize( int n, int nL){ return (((int)(n/nL))+1)*nL; }

    void cl_info(){
        printf( "sizeof(cl_mem    ) is %li bytes\n", sizeof(cl_mem) );
        printf( "sizeof(cl_kernel ) is %li bytes\n", sizeof(cl_kernel) );
        printf( "sizeof(cl_program) is %li bytes\n", sizeof(cl_program) );
        printf( "sizeof(cl_context) is %li bytes\n", sizeof(cl_context) );
        printf( "sizeof(cl_command_queue) is %li bytes\n", sizeof(cl_command_queue) );
        //printf( "sizeof(OCLtask   ) is %li bytes\n", sizeof(OCLtask) );
        //printf( "sizeof(OCLarg    ) is %li bytes\n", sizeof(OCLarg) );
        printf( "sizeof(OCLBuffer ) is %li bytes\n", sizeof(OCLBuffer) );
        //exit(0);
    }

    void print_device_params( cl_device_id device, bool bDetails=false ){
        cl_uint  ui;
        cl_ulong ul;
        size_t   sz;
        size_t   szs[4];
        const int nstrmax=1024;
        char     str[nstrmax];
        //clGetDeviceInfo(device, CL_DEVICE_NAME,                nstrmax,      str, NULL);   printf("DEVICE[%i,%i]: %s\n", i,j, str);
        //if(strstr(str, "NVIDIA") != NULL){ i_nvidia=i; }
        //clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR,              nstrmax,      str, NULL);   printf("\tVENDOR:      %s\n", str);
        printf("\t## VERSIONS: \n");
        clGetDeviceInfo(device, CL_DEVICE_VERSION,             nstrmax,      str, NULL);   printf("\tDEVICE_VERSION: %s\n", str);
        clGetDeviceInfo(device, CL_DRIVER_VERSION,             nstrmax,      str, NULL);   printf("\tDRIVER_VERSION: %s\n", str);
        clGetDeviceInfo(device, CL_DEVICE_OPENCL_C_VERSION,    nstrmax,      str, NULL);   printf("\tC_VERSION:      %s\n", str);
        printf("\t## PERFORMANCE: \n");
        clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS,        sizeof(uint),  &ui, NULL);  printf("\tMAX_COMPUTE_UNITS:    %d\n", ui  );
        clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY,      sizeof(uint),  &ui, NULL);  printf("\tMAX_CLOCK_FREQUENCY:  %d  MHz \n", ui   );
        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE,          sizeof(ulong), &ul, NULL);  printf("\tGLOBAL_MEM_SIZE:      %lu MB  \n", (ul/1024)/1024 );
        clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE,           sizeof(ulong), &ul, NULL);  printf("\tLOCAL_MEM_SIZE:       %lu kB  \n", ul/1024 );
        clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(ulong), &ul, NULL);  printf("\tCONSTANT_BUFFER_SIZE: %lu kB  \n", ul/1024 );
        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,    sizeof(ulong),  &ul, NULL); printf("\tGLOBAL_MEM_CACHE_SIZE:     %lu kB \n", (ul/1024) );
        clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE,sizeof(ulong),  &ul, NULL); printf("\tGLOBAL_MEM_CACHELINE_SIZE: %lu \n", ul );
        if(!bDetails) return;
        printf("\t## DETAILS: \n");
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,      sizeof(uint),    &ui, NULL); printf("\tMAX_WORK_ITEM_DIMENSIONS: %d  \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE,           sizeof(size_t),  &sz, NULL); printf("\tMAX_WORK_GROUP_SIZE:      %lu \n", sz );
        clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES,           sizeof(size_t)*4,szs, NULL); printf("\tMAX_WORK_ITEM_SIZES:      [%li,%li,%li] \n", szs[0],szs[1],szs[1] );
        clGetDeviceInfo(device, CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE,      sizeof(uint),    &ui, NULL); printf("\tMIN_DATA_TYPE_ALIGN_SIZE: %d  \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_IMAGE_MAX_BUFFER_SIZE,         sizeof(size_t),  &sz, NULL); printf("\tIMAGE_MAX_BUFFER_SIZE:    %d Mpix \n", sz/1024/1024 );
        clGetDeviceInfo(device, CL_DEVICE_IMAGE_MAX_ARRAY_SIZE,          sizeof(size_t),  &sz, NULL); printf("\tIMAGE_MAX_ARRAY_SIZE:     %d  \n", sz );
        clGetDeviceInfo(device, CL_DEVICE_QUEUE_ON_DEVICE_MAX_SIZE,      sizeof(uint),    &ui, NULL); printf("\tQUEUE_ON_DEVICE_MAX_SIZE:      %d   \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_QUEUE_ON_DEVICE_PREFERRED_SIZE,sizeof(uint),    &ui, NULL); printf("\tPREFERRED_QUEUE_SIZE:          %d \n", ui  );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR,   sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_CHAR:   %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT,  sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_SHORT:  %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT,    sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_INT:    %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG,   sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_LONG:   %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,  sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_FLOAT:  %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(uint),    &ui, NULL); printf("\tPREFERRED_VECTOR_WIDTH_DOUBLE: %d \n", ui );
        clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS,               nstrmax,       str, NULL);  printf("\tEXTENSIONS: %s\n", str);
        //printf("\n");
    }

    void printDeviceInfo( bool bDetails=false ){ print_device_params( device, bDetails ); }

    int print_devices( bool bDetails=false){
        int i, j;
        char*           value;
        size_t          valueSize;
        cl_uint         platformCount;
        cl_platform_id* platforms;
        cl_uint         deviceCount;
        cl_device_id*   devices;
        cl_uint         maxComputeUnits;
        // get all platforms
        clGetPlatformIDs(0, NULL, &platformCount);
        platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
        clGetPlatformIDs(platformCount, platforms, NULL);
        printf( "OpenCL platform count %i \n", platformCount );
        const int nstrmax=1024;
        char     str[nstrmax];
        cl_uint  ui;
        cl_ulong ul;
        size_t   sz;
        size_t   szs[4];
        int i_nvidia = -1;
        for (i = 0; i < platformCount; i++) {
            clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
            devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
            clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
            for (j = 0; j < deviceCount; j++) {
                clGetDeviceInfo(devices[j], CL_DEVICE_NAME,                nstrmax,      str, NULL);   printf("DEVICE[%i,%i]: %s\n", i,j, str);
                if(strstr(str, "NVIDIA") != NULL){ i_nvidia=i; }
                print_device_params( devices[j] );
                printf("\n");
            }
            free(devices);
        }
        free(platforms);
        return i_nvidia;
    }

    int copy( int iBufFrom, int iBufTo, int nbytes=-1, int src_offset=0, int dst_offset=0){
       if(nbytes<0){ nbytes=buffers[iBufFrom].byteSize(); int nbytes_=buffers[iBufTo].byteSize(); if(nbytes_<nbytes)nbytes=nbytes_;}
       //printf( "nbytes(-1) -> %i min(%i|%i) \n", nbytes, (int)buffers[iBufFrom].byteSize(), nbytes_ ); } 
       //printf( "OCLsystem::copy(%i[%i],%i[%i],n=%i)\n",iBufFrom,src_offset,iBufTo,dst_offset,nbytes  );
       int err = clEnqueueCopyBuffer( commands, buffers[iBufFrom].p_gpu, buffers[iBufTo].p_gpu, src_offset, dst_offset, nbytes, 0, 0, 0);
       OCL_checkError(err, "copy()");
       //finishRaw();
       return err;
    }

    int copyBuffToImage( int iBuff, int itex, size_t4 region, int src_offset=0 ){
        /*
        cl_int clEnqueueCopyBufferToImage(
            cl_command_queue command_queue,
            cl_mem src_buffer,
            cl_mem dst_image,
            size_t src_offset,
            const size_t* dst_origin,
            const size_t* region,
            cl_uint num_events_in_wait_list,
            const cl_event* event_wait_list,
            cl_event* event);
        */
        size_t offset[4]{0,0,0,0};
        //size_t region[4]{nx,ny,nz,0};
        printf( "copyBuffToImage() region(%li,%li,%li)\n", region.x, region.y, region.z  );
        int err = clEnqueueCopyBufferToImage( commands, buffers[iBuff].p_gpu, buffers[itex].p_gpu, src_offset, offset, (size_t*)&region, 0,0,0 );
        OCL_checkError(err, "copyBuffToImage()");
        return err;
    }

    /*
    int copyImageToBuffer( int iBuff, int itex, size_t4 region, int src_offset=0 ){
        printf( "copyBuffToImage() region(%li,%li,%li)\n", region.x, region.y, region.z  );
        err = clEnqueueCopyBufferToImage( commands, buffers[iBuff].p_gpu, buffers[itex].p_gpu, src_offset, offset, (size_t*)&region, 0,0,0 );
        OCL_checkError(err, "copyBuffToImage()");
        return err;
    }
    */

    void check_programSet (){ if(program ==0){ printf("ERROR OCLsystem program  not set \n"); exit(-1); } }
    void check_contextSet (){ if(context ==0){ printf("ERROR OCLsystem context  not set \n"); exit(-1); } }
    void check_deviceSet  (){ if(device  ==0){ printf("ERROR OCLsystem device   not set \n"); exit(-1); } }
    void check_commandsSet(){ if(commands==0){ printf("ERROR OCLsystem commands not set \n"); exit(-1); } }


    int autoDeviceChoice(int n, const cl_device_id* devices){
        int ipick = -1;
        for(int i=0; i<n; i++){
            char name[MAX_INFO_STRING];
            getDeviceName(devices[i], name);
            printf("Device[%i]: %s\n", i, name);
            if(strstr(name, "NVIDIA") != NULL){ ipick=i; }
        }
        if(ipick<0){ printf( "OCLsystem::autoDeviceChoice() WARRNING nVidia not found!\n" ); ipick=0; }
        return ipick;
    }

    int initOCL(int deviceIndex=-1 ){
        int err=0;
        //cl_info(); exit(0);
        //cl_uint deviceIndex = 0;
        //parseArguments(argc, argv, &deviceIndex);
        cl_device_id devices[MAX_DEVICES];
        unsigned numDevices = getDeviceList(devices);
        if( deviceIndex<0 ){
            deviceIndex = autoDeviceChoice(numDevices, devices );
        }
        if ((deviceIndex >= numDevices)||(deviceIndex<0)){  printf("ERROR in OCLsystem::initOCL() deviceIndex(%i)>=numDevices(%i)\n", deviceIndex, numDevices ); exit(0); }
        device = devices[deviceIndex];
        char name[MAX_INFO_STRING];
        getDeviceName(device, name);
        printf("\nUsing OpenCL device: %s\n", name);
        context  = clCreateContext(0, 1, &device, NULL, NULL, &err);  OCL_checkError(err, "Creating context");
        //commands = clCreateCommandQueue(context, device, 0, &err);    OCL_checkError(err, "Creating command queue");   // DEPRECATED
        //cl_uint maxQueueSize = 450000;
        //cl_queue_properties prop[] = { CL_QUEUE_PROPERTIES, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_ON_DEVICE | CL_QUEUE_ON_DEVICE_DEFAULT,  CL_QUEUE_SIZE, maxQueueSize, 0 };
        cl_command_queue_properties prop = 0;
        //prop |= CL_QUEUE_PROFILING_ENABLE;
        //prop |= CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
        //commands = clCreateCommandQueueWithProperties( context, device, &prop, &err );
        commands = clCreateCommandQueueWithProperties( context, device, NULL, &err );   // if properties are NULL defaults are used 
        OCL_checkError(err, "Creating command queue");
        return err;
    }
    
    int newTask( const char * name, cl_program program_=0, size_t dim=1, size_t4 global=size_t4{0,0,0,0}, size_t4 local={1,1,1,1} ){
        if(program_==0){ 
            check_programSet();
            program_=program;
        }
        int iker = newKernel( name, program_ );
        OCLtask* task = new OCLtask( this, iker, dim, global, local );
        tasks.push_back( task ); 
        int i = tasks.size()-1;
        task_dict.insert({name,i});
        kernel_dict.insert( { name, i } );
        return i;
    }

    int newKernel( const char * name, cl_program program_=0 ){
        if(program_==0){ 
            check_programSet();
            program_=program;
        }
        //printf( "newKernel() program %li %li \n", (long)program_, (long)program );
        int err=0; kernels.push_back( clCreateKernel( program_, name, &err ) );  OCL_checkError(err, "newKernel"); 
        int i = kernels.size()-1;
        kernel_dict.insert( { name, i } );
        return i;
    }

    int newBuffer( const char* name, size_t n, size_t typesize, void* p_cpu=0, cl_mem_flags flags=CL_MEM_READ_WRITE ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, n, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffer_dict.insert( { name, i } );
        int err=buffers[i].initOnGPU(context); OCL_checkError__(err, "newBuffer",i,name);
        return i;
    }

    int newBufferImage2D( const char* name, size_t nx, size_t ny, size_t typesize, void * p_cpu, cl_mem_flags flags, cl_image_format imageFormat ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, nx*ny, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffer_dict.insert( { name, i } );
        buffers[i].img_dims    = 2;
        buffers[i].nImg[0]     = nx;
        buffers[i].nImg[1]     = ny;
        buffers[i].imageFormat = imageFormat;
        int err=buffers[i].initOnGPUImage(context); OCL_checkError__(err, "newBufferImage2D",i,name);
        return i;
    }

    int newBufferImage3D( const char* name, size_t nx, size_t ny, size_t nz, size_t typesize, void * p_cpu, cl_mem_flags flags, cl_image_format imageFormat ){
        check_contextSet();
        buffers.push_back( OCLBuffer( name, nx*ny, typesize, p_cpu, flags ) );
        int i=buffers.size()-1;
        buffer_dict.insert( { name, i } );
        buffers[i].img_dims    = 3;
        buffers[i].nImg[0]     = nx;
        buffers[i].nImg[1]     = ny;
        buffers[i].nImg[2]     = nz;
        buffers[i].imageFormat = imageFormat;
        //printf( "DEBUG newBufferImage2D() flags %li %li \n", flags, flags );
        //printf( "newBufferImage3D buffers[%i].img_dims %i \n", i, buffers[i].img_dims  );
        int err=buffers[i].initOnGPUImage(context); OCL_checkError__(err, "newBufferImage3D",i,name );
        return i;
    }

    int initBuffers   (){ int err = CL_SUCCESS; for(size_t i=0; i<buffers.size(); i++){  err |= buffers[i].initOnGPU ( context );     }; return err; }
    //int releaseBuffers(){ for(int i=0; i<buffers; i++){ clReleaseMemObject(buffers[i].p_gpu); } }

    char * getKernelSource(const char *filename){
        FILE *file = fopen(filename, "r");
        if (!file){ fprintf(stderr, "Error: Could not open kernel source file\n"); exit(-1); }
        fseek(file, 0, SEEK_END);
        int len = ftell(file) + 1;
        rewind(file);
        char *source = (char *)calloc(sizeof(char), len);
        //char* source = new char[len];
        if (!source){ fprintf(stderr, "Error: Could not allocate memory for source string\n"); exit(-1); }
        fread(source, sizeof(char), len, file);
        fclose(file);
        return source;
    }

    int buildProgram( const char * fname, cl_program& program_ ){       // TODO : newProgram instead ?
        int err=0;
        char * kernelsource = getKernelSource( fname );
        program_ = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
        char tmpstr[1024];
        sprintf(tmpstr,"Creating program with %s", fname);
        OCL_checkError(err, tmpstr);
        free(kernelsource);
        err =      clBuildProgram(program_, 0,         NULL,      "-I. -cl-std=CL2.0", NULL, NULL);
        //free(kernelsource);     // Why it crashes ?
        if (err != CL_SUCCESS){
            printf( " ERROR in clBuildProgram %s \n", fname);
            OCL_buildProgramFailure( program_, device );
            exit(0);
            return -1;
        }
        //delete [] kernelsource; // TODO ??????
        //free(kernelsource);     // Why it crashes ?
        return err;
    }
    int buildProgram( const char * fname ){ return buildProgram( fname, program ); }
 
    inline int upload  (int i, const void* cpu_data, int n=-1,int i0=0 ){ return buffers[i].toGPU  (commands,cpu_data,n,i0); };
    inline int download(int i,       void* cpu_data, int n=-1,int i0=0 ){ return buffers[i].fromGPU(commands,cpu_data,n,i0); };
    //inline int upload  (int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].toGPU  (commands); };
    //inline int download(int i, void* p_cpu ){ buffers[i].p_cpu=p_cpu; return buffers[i].fromGPU(commands); };

    inline int upload  (int i){ return buffers[i].toGPU  (commands); };
    inline int download(int i){ return buffers[i].fromGPU(commands); };

    //inline int copy_raw    (int from, int to, int from0, int to0, int n){ return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,from0,to0,n,0,NULL,NULL); };
    //inline int copyBuff(int from, int to                           ){ int n=buffers[from].n; int n_=buffers[to].n; if(n_<n)n=n_; return clEnqueueCopyBuffer(commands,buffers[from].p_gpu,buffers[to].p_gpu,0,0,n,0,NULL,NULL); };

    void useKernel( int ikernel){ current_kernel=kernels[ikernel]; argCounter=0; };
    void useTask( OCLtask* task ){ useKernel( task->ikernel );  currentTask=task; };

    // void errArg(cl_int err){ currentTask 
    //     if (err != CL_SUCCESS){ 
    //         OCL_checkError_(err,"useArg",i,currentTask->n); 
    //     }
    // }

    //int  useArg( cl_mem ibuff,              int i=-1 ){ if(i<0){i=argCounter;argCounter++;} int err=clSetKernelArg( current_kernel, i, sizeof(cl_mem), &(ibuff) ); printf("useArg[%i]\n",i); OCL_checkError_(err,"useArg",i); return err; };
    int  useArg( int    i_arg,              int i=-1 ){ if(i<0){i=argCounter;argCounter++;} int err=clSetKernelArg( current_kernel, i, sizeof(int),    &(i_arg) );                OCL_error_warn_(err,"useArg",i);     return err; };
    int  useArg( float  f_arg,              int i=-1 ){ if(i<0){i=argCounter;argCounter++;} int err=clSetKernelArg( current_kernel, i, sizeof(float),  &(f_arg) );                OCL_error_warn_(err,"useArg",i);     return err; };
    int  useArg_( void*  buff , int nbytes, int i=-1 ){ if(i<0){i=argCounter;argCounter++;} int err=clSetKernelArg( current_kernel, i, nbytes,           buff   );                OCL_error_warn_(err,"useArg_",i);    return err; };
    int  useArgBuff( int ibuff,             int i=-1 ){ if(i<0){i=argCounter;argCounter++;} int err=clSetKernelArg( current_kernel, i, sizeof(cl_mem), &(buffers[ibuff].p_gpu) ); OCL_error_warn_(err,"useArgBuff",i); return err; };
    int enque( size_t dim, const size_t* global, const size_t* local, int ikernel=-1 ){ 
        cl_kernel kernel;
        if(ikernel<0){kernel=current_kernel;}else{ kernel = kernels[ikernel]; }; 
        //printf( "OCLsystem::enque() dim %li global[%li,%li,%li] local[%li,%li,%li] \n", dim, global[0],global[1],global[2], local[0],local[1],local[2] );
        return clEnqueueNDRangeKernel( commands, kernel, dim, NULL, global, local, 0, NULL, NULL ); 
    }

    int enque( size_t dim, size_t4 global, size_t4 local, int ikernel=-1 ){  return enque(dim, (size_t*)&global, (size_t*)&local, ikernel); }
    int enque( size_t dim, size_t4 global, int ikernel=-1 ){  return enque(dim, (size_t*)&global, NULL, ikernel); }

    int enqueTask( int i ){ return tasks[i]->enque(); };
    int enqueTask( const char* name ){ return enqueTask( task_dict[name] ); }

    int download(){
        int err = CL_SUCCESS;
        for(size_t i=0; i<buffers.size(); i++ ){
            if( buffers[i].read_on_finish ){
                //printf("finish : reading buff %i \n", i);
                err |= buffers[i].fromGPU( commands );
            }
        }
        return err;
    }

    int finishRaw(){
        int err = clFinish(commands); return err;
    }

    int finish(){
        int err=0;
        err = clFinish(commands);
        err |= download();
        return err;
    }

    void release_OCL( cl_program program=0 ){
        printf( "OCL_DEF::release_OCL() --- START \n" );
        if(program==0) program=this->program;
        printf( "OCL_DEF::release_OCL() kernels \n" );
        for(size_t i=0; i<kernels.size(); i++){ clReleaseKernel   (kernels[i]      ); }
        printf( "OCL_DEF::release_OCL() buffers \n" );
        for(size_t i=0; i<buffers.size(); i++){ clReleaseMemObject(buffers[i].p_gpu); }
        printf( "OCL_DEF::release_OCL() program \n" );
        clReleaseProgram(program);   
        printf( "OCL_DEF::release_OCL() commands \n" );
        //clReleaseKernel(kernel);
        clReleaseCommandQueue(commands);
        printf( "OCL_DEF::release_OCL() context \n" );
        clReleaseContext(context);
        printf( "OCL_DEF::release_OCL() --- DONE \n" );
    }
    ~OCLsystem(){ 
        release_OCL();
    }

    void printBuffers(){
        printf("OCLsystem::printBuffers()\n");
        for(int i=0; i<buffers.size(); i++ ){
            const OCLBuffer& b = buffers[i];
            printf( "buffers[%i](%s) %li*%li img_dims %i\n", i, b.name.c_str(), b.n, b.typesize, b.img_dims );
        } 
    }
}; // class OCLsystem

//=======================================================================
//=======================================================================

int OCLtask::enque_raw(  ){
    //printf("enque_raw ikernel %li idim %li global(%li,%li,%li) local(%li,%li,%li)\n", ikernel, dim, global.x,global.y,global.z, local.x,local.y,local.z );    
    if(local.x==0){ return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, (size_t*)&global, NULL,            0, NULL, NULL );   }
    else{           return clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, (size_t*)&global, (size_t*)&local, 0, NULL, NULL );   }
}

int OCLtask::useArgs( ){
    int err = CL_SUCCESS;
    cl_kernel kernel = cl->kernels[ikernel];
    for(size_t i=0; i<args.size(); i++){
        OCLarg& arg = args[i];
        switch(arg.kind){
            case OCL_BUFF:
                //printf( "DEBUG buffArg args[%li] ibuff %i p_gpu %li '%s' \n", i, arg.i, (long)cl->buffers[arg.i].p_gpu, cl->buffers[arg.i].name );
                err |= clSetKernelArg( kernel, i, sizeof(cl_mem), &(cl->buffers[arg.i].p_gpu) );             break;
            //case OCL_BUFF:  err |= cl->buffers[arg.i].setAsArg( kernel, i );                               break;
            case OCL_INT:   err |= clSetKernelArg( kernel, i, sizeof(int)  , &(arg.i) );                     break;
            case OCL_FLOAT: err |= clSetKernelArg( kernel, i, sizeof(float), &(arg.f) );                     break;
            case OCL_LBUFF: err |= clSetKernelArg( kernel, i, arg.i,          NULL    );                     break;
            case OCL_PTR:   err |= clSetKernelArg( kernel, i, arg.nbytes,     arg.ptr );                     break;
        }
        //OCL_checkError(err, "setAsArg", i );
        if(bOCLCheckError)OCL_check_error(err,"setAsArg",__FILE__,__LINE__,i);
    }
    return err;
}

void OCLtask::print_arg_list( ){
    printf("DEBUG print_arg_list \n" );
    printf("kernel[narg=%li]( \n", args.size() );
    for(size_t i=0; i<args.size(); i++){
        printf( "arg %li \n", i );
        switch(args[i].kind){
            case OCL_INT:   printf( "[%li]int %i, ",     i, args[i].i ); break;
            case OCL_FLOAT: printf( "[%li]float %g, ",   i, args[i].f ); break;
            //case OCL_BUFF:  printf( "[%li]buff[%i]:%s, ",i, task.args[i].i, buffers[task.args[i].i].name.c_str() );   break;
            case OCL_BUFF:  printf( "[%li]buff[%i] ",i, args[i].i );   break;
            default:        printf( "[%li]arg(type=%i,val=%i) ",  i, args[i].kind, args[i].i );   break;
        }
    }
    printf(")\n");
}

#endif
