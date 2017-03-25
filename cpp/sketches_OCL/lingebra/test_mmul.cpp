#include <cstdio>
#include <math.h>

#include "fastmath.h"
#include "testUtils.h"
//#include "Lingebra.h"

#include "OCLerrors.h"
#include "OCL_device_picker.h"
#include "OCL.h"

// ==================== Matrix and vector utilities

float error( int n, float * xs, float * x0s ){
    float sum = 0;
    for(int i=0; i<n; i++){ float d = xs[i]-x0s[i]; sum+=d*d; }
    return sum;
}

float setArrayRandom( int n, float xmin, float xmax, float * out ){
    float xrange = xmax - xmin;
    for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	}
}

float setArray( int n, float val, float * out ){  for (int i=0; i<n; i++ ){ out[i] = val; } }

void mmul_cpu(int N, float *A, float *B, float *C){
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float tmp = 0.0f;
            for (int k = 0; k < N; k++) {  tmp += A[i*N+k] * B[k*N+j];  }
            C[i*N+j] = tmp;
        }
    }
}

void mat2file(FILE *f, int n, float* A){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(f," %10.10e", A[i*n+j] );
        }
        fprintf(f,"\n");
    }
}

// ==================== Globals

OCLsystem cl;

cl_int err;
int nrep = 1;
//int n    = 32;
//int n    = 64;
//int n    = 128;
//int n    = 256;
int n    = 1024;
int n2   = n*n;

float * A    = new float[ n2 ];
float * B    = new float[ n2 ];
float * C    = new float[ n2 ];
float * Cref = new float[ n2 ];

OCLtask * task1;

long t1;
double fticks1;
double fticks2;

// ==================== Task derived types

class Task_Mmul_row : public OCLtask {
    public:
    int n;
    virtual int enque( ){
        int err;
        cl_kernel kernel = cl->kernels[ikernel];
        err = clSetKernelArg( kernel, 0, sizeof(int),  &n );
        err = cl->buffers[0].setAsArg( kernel, 1 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[1].setAsArg( kernel, 2 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[2].setAsArg( kernel, 3 );     OCL_checkError(err, "setAsArg");
        //err = clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL ); OCL_checkError(err, "clEnqueueNDRangeKernel");
        err = enque_raw( );                             OCL_checkError(err, "enque_raw");
    }
};

class Task_Mmul_privB : public OCLtask {
    public:
    int n;
    virtual int enque( ){
        int err;
        cl_kernel kernel = cl->kernels[ikernel];
        err = clSetKernelArg( kernel, 0, sizeof(int),  &n );
        err = cl->buffers[0].setAsArg( kernel, 1 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[1].setAsArg( kernel, 2 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[2].setAsArg( kernel, 3 );     OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 4, sizeof(float) * n,  NULL );  OCL_checkError(err, "setAsArg");
        //err = clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL ); OCL_checkError(err, "clEnqueueNDRangeKernel");
        err = enque_raw( );                             OCL_checkError(err, "enque_raw");
    }
};

class Task_Mmul_local : public OCLtask {
    public:
    int n;
    virtual int enque( ){
        int err;
        cl_kernel kernel = cl->kernels[ikernel];
        err = cl->buffers[0].setAsArg( kernel, 0 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[1].setAsArg( kernel, 1 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[2].setAsArg( kernel, 2 );     OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 3, sizeof(int),  &n );   OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 4, sizeof(int),  &n );   OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 5, sizeof(int),  &n );   OCL_checkError(err, "setAsArg");
        //err = clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL ); OCL_checkError(err, "clEnqueueNDRangeKernel");
        err = enque_raw( );                             OCL_checkError(err, "enque_raw");
    }
};

class Task_Mmul_block : public OCLtask {
    public:
    int n;
    virtual int enque( ){
        int err;
        int block_size = 16;
        cl_kernel kernel = cl->kernels[ikernel];
        err = clSetKernelArg( kernel, 0, sizeof(int),  &n );
        err = cl->buffers[0].setAsArg( kernel, 1 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[1].setAsArg( kernel, 2 );     OCL_checkError(err, "setAsArg");
        err = cl->buffers[2].setAsArg( kernel, 3 );     OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 4, sizeof(float) * block_size * block_size,  NULL );  OCL_checkError(err, "setAsArg");
        err = clSetKernelArg( kernel, 5, sizeof(float) * block_size * block_size,  NULL );  OCL_checkError(err, "setAsArg");
        printf( "global_sz (%i,%i) (%i,%i) \n", global[0], global[1], local[0], local[1] );
        err = clEnqueueNDRangeKernel( cl->commands, cl->kernels[ikernel], dim, NULL, global, local, 0, NULL, NULL ); OCL_checkError(err, "clEnqueueNDRangeKernel");
        //err = enque_raw( );                             OCL_checkError(err, "enque_raw");
    }
};

// ==================== Main

void runCL_direct(){
    printf("DEBUG = 1\n");
    cl_mem A_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, A,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem B_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, B,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem C_clmem = clCreateBuffer( cl.context, CL_MEM_WRITE_ONLY                       , sizeof(float) * n2, NULL  , &err );   OCL_checkError(err, "clCreateBuffer");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmuls.cl" );                          OCL_checkError(err, "cl.buildProgram");
    cl_kernel kernel = clCreateKernel( cl.program, "mmul_row", &err );    OCL_checkError(err, "clCreateKernel");
    //cl_kernel kernel = clCreateKernel( cl.program, "mmul_block", &err );    OCL_checkError(err, "clCreateKernel");
    printf("DEBUG = 4\n");
    size_t global_sz[2] = {n2, n }; size_t* local_sz    = NULL;   //  mmul_row
    //size_t global_sz[2] = {n, n }; size_t  local_sz[2] = {16,16};   //  mmul_block
    t1 = getCPUticks();
    for (int irep = 0; irep < nrep; irep++){
        err =  clSetKernelArg( kernel, 0, sizeof(int),    &n       );        OCL_checkError(err, "clSetKernelArg");
        err |= clSetKernelArg( kernel, 1, sizeof(cl_mem), &A_clmem );        OCL_checkError(err, "clSetKernelArg");
        err |= clSetKernelArg( kernel, 2, sizeof(cl_mem), &B_clmem );        OCL_checkError(err, "clSetKernelArg");
        err |= clSetKernelArg( kernel, 3, sizeof(cl_mem), &C_clmem );        OCL_checkError(err, "clSetKernelArg");
        //err |= clSetKernelArg( kernel, 4, sizeof(float) * 16 * 16,  NULL );  OCL_checkError(err, "clSetKernelArg");
        //err |= clSetKernelArg( kernel, 5, sizeof(float) * 16 * 16,  NULL );  OCL_checkError(err, "clSetKernelArg");
        printf("DEBUG 4 (%i)\n", irep );
        err = clEnqueueNDRangeKernel( cl.commands, kernel, 1, NULL, global_sz, local_sz, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel");  // mmul_row
        //err = clEnqueueNDRangeKernel( cl.commands, kernel, 2, NULL, global_sz, local_sz, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel"); // mmul_block
    }
    printf("DEBUG = 5\n");
    err = clFinish           ( cl.commands );                                                              OCL_checkError(err, "clFinish");
    fticks1 = (getCPUticks()-t1);
    printf("DEBUG = 6\n");
    err = clEnqueueReadBuffer( cl.commands, C_clmem, CL_TRUE, 0, sizeof(float) * n2, C, 0, NULL, NULL);    OCL_checkError(err, "clEnqueueReadBuffer");
    fticks2 = (getCPUticks()-t1);
}

void runCL_auto(){
    printf("DEBUG 2\n");
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), A, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), B, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), C, CL_MEM_WRITE_ONLY ) );
    cl.buffers[2].read_on_finish = true;
    printf("DEBUG 2\n");
    err = cl.initBuffers();                                                    OCL_checkError(err, "clinitBuffers");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmuls.cl" );                               OCL_checkError(err, "cl.buildProgram");
    // ---- mmul_row
    //cl.kernels.push_back( clCreateKernel( cl.program, "mmul_row", &err ) );    OCL_checkError(err, "clCreateKernel");
    //Task_Mmul_row   * task_ = new Task_Mmul_row();   task_->setup( &cl, cl.kernels.size()-1, 1, n2, 0    );
    // ---- mmul_row_priv
    //cl.kernels.push_back( clCreateKernel( cl.program, "mmul_row_priv", &err ) );    OCL_checkError(err, "clCreateKernel");
    //Task_Mmul_row * task_ = new Task_Mmul_row();   task_->setup( &cl, cl.kernels.size()-1, 1, n2, 0    );
    // ---- mmul_row_priv_block
    //cl.kernels.push_back( clCreateKernel( cl.program, "mmul_row_priv_block", &err ) );    OCL_checkError(err, "clCreateKernel");
    //Task_Mmul_privB   * task_ = new Task_Mmul_privB();   task_->setup( &cl, cl.kernels.size()-1, 1, n2, n/16    );
    // ---- mmul_local
    cl.kernels.push_back( clCreateKernel( cl.program, "mmul_local", &err ) );    OCL_checkError(err, "clCreateKernel");
    Task_Mmul_local   * task_ = new Task_Mmul_local();   task_->setup( &cl, cl.kernels.size()-1, 2, n, 16    );
    // ---- mmul_block
    //cl.kernels.push_back( clCreateKernel( cl.program, "mmul_block", &err ) );    OCL_checkError(err, "clCreateKernel");
    //Task_Mmul_block * task_ = new Task_Mmul_block(); task_->setup( &cl, cl.kernels.size()-1, 2, n, 16 );
    task_->n = n;
    task1    = task_;
    printf("DEBUG = 4\n");
    t1 = getCPUticks();
    for (int irep = 0; irep < nrep; irep++){
        task1->enque();
    }
    printf("DEBUG = 5\n");
    //err = cl.finish();            OCL_checkError(err, "cl.finish");
    err = clFinish(cl.commands);    OCL_checkError(err, "clFinish");
    fticks1 = (getCPUticks()-t1);
    err = cl.download();            OCL_checkError(err, "cl.download");
    fticks2 = (getCPUticks()-t1);
}

int main(){
    //printf("DEBUG 0.1\n");

    setArrayRandom( n2, -1.0, 1.0, A );       //printf("DEBUG 0.1\n");
    setArrayRandom( n2, -1.0, 1.0, B );       //printf("DEBUG 0.2\n");

    printf(" matrix size = %i \n", n );
    printf("DEBUG 1\n");


    cl.init();

    //runCL_direct();
    runCL_auto();
    printf( "TIME( GPU          ) %g [ticks] %g [tick/op]\n", fticks1, fticks1/(n*n*n) );
    printf( "TIME( GPU+download ) %g [ticks] %g [tick/op]\n", fticks2, fticks2/(n*n*n) );

    printf( "==== check on CPU \n" );

    setArray      ( n2,  0.0,      Cref );    //printf("DEBUG 0.3\n");
    //setArray      ( cl.buffers[2].n,  0.0,       (float*)cl.buffers[2].C_cpu );
    t1 = getCPUticks();
    mmul_cpu( n, A, B, Cref );                // printf("DEBUG 0.4\n");
    double fticks_cpu = (getCPUticks()-t1);
    printf( "TIME( CPU          ) %g [ticks] %g [tick/op]\n", fticks_cpu, fticks_cpu/(n*n*n) );
    printf( "time CPU/GPU %g +download %g \n", fticks_cpu/fticks1, fticks_cpu/fticks2 );

    float delta = error( n2, C, Cref);
    printf(  "numerical error = %g \n", delta  );

    printf( "==== save matrices for control ... " );
    FILE *fout;
    fout=fopen("A.dat",    "w"); mat2file( fout, n, A   ); fclose(fout);
    fout=fopen("B.dat",    "w"); mat2file( fout, n, B   ); fclose(fout);
    fout=fopen("C.dat",    "w"); mat2file( fout, n, C   ); fclose(fout);
    fout=fopen("Cref.dat", "w"); mat2file( fout, n, Cref); fclose(fout);
    printf( "DONE\n" );

}
