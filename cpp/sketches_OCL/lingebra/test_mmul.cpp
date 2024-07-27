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

void setArrayRandom( int n, float xmin, float xmax, float * out ){
    float xrange = xmax - xmin;
    for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	}
}

void setArray( int n, float val, float * out ){  for (int i=0; i<n; i++ ){ out[i] = val; } }

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

// ==================== Main

void runCL_direct(){
    printf("DEBUG = 1\n");
    cl_mem A_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, A,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem B_clmem = clCreateBuffer( cl.context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR , sizeof(float) * n2, B,      &err );   OCL_checkError(err, "clCreateBuffer");
    cl_mem C_clmem = clCreateBuffer( cl.context, CL_MEM_WRITE_ONLY                       , sizeof(float) * n2, NULL  , &err );   OCL_checkError(err, "clCreateBuffer");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmuls.cl" );                                   OCL_checkError(err, "cl.buildProgram");
    cl_kernel kernel = clCreateKernel( cl.program, "mmul_row", &err );        OCL_checkError(err, "clCreateKernel");
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
    cl.newBuffer( "A", n2, sizeof(float), A, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "B", n2, sizeof(float), B, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR );
    cl.newBuffer( "C", n2, sizeof(float), C, CL_MEM_WRITE_ONLY );
    cl.buffers[2].read_on_finish = true;
    //OCL_checkError(err, "clinitBuffers");
    printf("DEBUG = 3\n");
    err = cl.buildProgram( "cl/mmuls.cl" );                               OCL_checkError(err, "cl.buildProgram");
    //int iker = cl.newKernel("mmul_row");            task1 = new OCLtask( &cl, iker, 1, n2, 0    ); task1->args = { INTarg(n), BUFFarg(0), BUFFarg(1), BUFFarg(2) };
    //int iker = cl.newKernel("mmul_row_priv");       task1 = new OCLtask( &cl, iker, 1, n2, 0    ); task1->args = { INTarg(n), BUFFarg(0), BUFFarg(1), BUFFarg(2) };
    //int iker = cl.newKernel("mmul_row_priv_block"); task1 = new OCLtask( &cl, iker, 1, n2, n/16 ); task1->args = { INTarg(n), BUFFarg(0), BUFFarg(1), BUFFarg(2), LBUFFarg(n*sizeof(float)) };
    int iker = cl.newKernel("mmul_local");            task1 = new OCLtask( &cl, iker, 2, n, 16    ); task1->args = {            BUFFarg(0), BUFFarg(1), BUFFarg(2), INTarg(n),INTarg(n),INTarg(n) };
    //int iker = cl.newKernel("mmul_block");          task1 = new OCLtask( &cl, iker, 2, n, 16    ); task1->args = { INTarg(n), BUFFarg(0), BUFFarg(1), BUFFarg(2), LBUFFarg(16*16*sizeof(float)), LBUFFarg(16*16*sizeof(float)) };

    task1->print_arg_list();

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


    cl.initOCL();

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
