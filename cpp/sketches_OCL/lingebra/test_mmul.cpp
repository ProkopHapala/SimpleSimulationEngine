#include <cstdio>
#include <math.h>

#include "fastmath.h"
//#include "Lingebra.h"

#include "OCLerrors.h"
#include "OCL_device_picker.h"
#include "OCL.h"

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
int n    = 32;
int n2   = n*n;

float * A    = new float[ n2 ];
float * B    = new float[ n2 ];
float * C    = new float[ n2 ];
float * Cref = new float[ n2 ];

// ==================== Main

void runCL_direct(){
    printf("DEBUG = 1\n");
    cl_mem A_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, A,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem B_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, B,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem C_clmem = clCreateBuffer( cl.context, CL_MEM_WRITE_ONLY                       , sizeof(float) * n2, NULL  , &err );   OCL_checkError(err, "clCreateBuffer");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmul/C_row.cl" );
    cl_kernel kernel = clCreateKernel( cl.program, "mmul", &err );          OCL_checkError(err, "clCreateKernel");
    printf("DEBUG = 4\n");
    size_t global_sz[2] = {n2, n };
    //size_t local_sz[]   = {};
    size_t* local_sz    = NULL;
    for (int irep = 0; irep < nrep; irep++){
        err =  clSetKernelArg( kernel, 0, sizeof(int),    &n       );
        err |= clSetKernelArg( kernel, 1, sizeof(cl_mem), &A_clmem );
        err |= clSetKernelArg( kernel, 2, sizeof(cl_mem), &B_clmem );
        err |= clSetKernelArg( kernel, 3, sizeof(cl_mem), &B_clmem );
        OCL_checkError(err, "clSetKernelArg");
        printf("DEBUG 4 (%i)\n", irep );
        //err = clEnqueueNDRangeKernel( cl.commands, cl.kernels[0].kernel, 1, NULL, cl.kernels[0].global, cl.kernels[0].local, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel");
        err = clEnqueueNDRangeKernel( cl.commands, kernel, 1, NULL, global_sz, local_sz, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel");
    }
    printf("DEBUG = 5\n");
    err = clFinish           ( cl.commands );                                                                                            OCL_checkError(err, "clFinish");
    printf("DEBUG = 6\n");
    err = clEnqueueReadBuffer( cl.commands, C_clmem, CL_TRUE, 0, sizeof(float) * n2, C, 0, NULL, NULL);    OCL_checkError(err, "clEnqueueReadBuffer");
}

void runCL_semi_direct(){
    printf("DEBUG = 1\n");
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), A, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), B, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), C, CL_MEM_WRITE_ONLY ) );
    cl.buffers[2].read_on_finish = true;
    printf("DEBUG = 2\n");
    cl.buffers[0].p_gpu = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * n2, cl.buffers[0].p_cpu, &err ); OCL_checkError(err, "clCreateBuffer");
    cl.buffers[1].p_gpu = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * n2, cl.buffers[1].p_cpu, &err ); OCL_checkError(err, "clCreateBuffer");
    cl.buffers[2].p_gpu = clCreateBuffer( cl.context, CL_MEM_WRITE_ONLY                      , sizeof(float) * n2, NULL               , &err ); OCL_checkError(err, "clCreateBuffer");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmul/C_row.cl" );                                                  OCL_checkError(err, "buildProgram");
    cl.kernels.push_back( OCLKernel( clCreateKernel( cl.program, "mmul", &err ), 1, n2, 0 ) );    OCL_checkError(err, "clCreateKernel");
    printf("DEBUG = 4\n");
    for (int irep = 0; irep < nrep; irep++){
        err =  clSetKernelArg( cl.kernels[0].kernel, 0, sizeof(int),    &n                   );
        err |= clSetKernelArg( cl.kernels[0].kernel, 1, sizeof(cl_mem), &cl.buffers[0].p_gpu );
        err |= clSetKernelArg( cl.kernels[0].kernel, 2, sizeof(cl_mem), &cl.buffers[1].p_gpu );
        err |= clSetKernelArg( cl.kernels[0].kernel, 3, sizeof(cl_mem), &cl.buffers[2].p_gpu );
        OCL_checkError(err, "clSetKernelArg");
        printf("DEBUG 4 (%i)\n", irep );
        //err = clEnqueueNDRangeKernel( cl.commands, cl.kernels[0].kernel, 1, NULL, cl.kernels[0].global, cl.kernels[0].local, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel");
        err = clEnqueueNDRangeKernel( cl.commands, cl.kernels[0].kernel, 1, NULL, cl.kernels[0].global, NULL, 0, NULL, NULL);   OCL_checkError(err, "clEnqueueNDRangeKernel");
    }
    printf("DEBUG = 5\n");
    err = clFinish              ( cl.commands );                                                                                           OCL_checkError(err, "clFinish");
    printf("DEBUG = 6\n");
    err = clEnqueueReadBuffer   ( cl.commands, cl.buffers[2].p_gpu, CL_TRUE, 0, sizeof(float) * n2, cl.buffers[2].p_cpu, 0, NULL, NULL);    OCL_checkError(err, "clEnqueueReadBuffer");
}

void runCL_auto(){
    printf("DEBUG 2\n");
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), A, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), B, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR ) );
    cl.buffers.push_back( OCLBuffer( n2, sizeof(float), C, CL_MEM_WRITE_ONLY ) );
    cl.buffers[2].read_on_finish = true;
    printf("DEBUG 2\n");
    err = cl.initBuffers();                                                                       OCL_checkError(err, "clinitBuffers");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmul/C_row.cl" );                                                  OCL_checkError(err, "cl.buildProgram");
    cl.kernels.push_back( OCLKernel( clCreateKernel( cl.program, "mmul", &err ), 1, n2, 0 ) );    OCL_checkError(err, "clCreateKernel");
    printf("DEBUG = 4\n");
    for (int irep = 0; irep < nrep; irep++){
        err  = clSetKernelArg( cl.kernels[0].kernel, 0, sizeof(int), &n );
        err |= cl.buffers[0].setAsArg( cl.kernels[0].kernel, 1 );    OCL_checkError(err, "clSetKernelArg");
        err |= cl.buffers[1].setAsArg( cl.kernels[0].kernel, 2 );    OCL_checkError(err, "clSetKernelArg");
        err |= cl.buffers[2].setAsArg( cl.kernels[0].kernel, 3 );    OCL_checkError(err, "clSetKernelArg");
        printf("DEBUG 4 (%i)\n", irep );
        cl.kernels[0].enque( cl.commands );
    }
    printf("DEBUG = 5\n");
    err = cl.finish();                  OCL_checkError(err, "cl.finish");
}

int main(){

    //printf("DEBUG 0.1\n");

    setArrayRandom( n2, -1.0, 1.0, A );       //printf("DEBUG 0.1\n");
    setArrayRandom( n2, -1.0, 1.0, B );       //printf("DEBUG 0.2\n");
    setArray      ( n2,  0.0,      Cref );    //printf("DEBUG 0.3\n");
    //setArray      ( cl.buffers[2].n,  0.0,       (float*)cl.buffers[2].C_cpu );
    mmul_cpu( n, A, B, Cref );  printf("DEBUG 0.4\n");

    FILE *fout;
    fout=fopen("A.dat",    "w"); mat2file( fout, n, A   ); fclose(fout);
    fout=fopen("B.dat",    "w"); mat2file( fout, n, B   ); fclose(fout);
    fout=fopen("Cref.dat", "w"); mat2file( fout, n, Cref); fclose(fout);

    printf("DEBUG 1\n");

    cl.init();

    //runCL_direct();
    //runCL_semi_direct();
    runCL_auto();


    fout=fopen("C.dat", "w"); mat2file( fout, n, C); fclose(fout);
    float delta = error( n2, C, Cref);

    printf(  "numerical error = %g \n", delta  );

}
