#include <cstdio>
#include <math.h>

#include "fastmath.h"
//#include "Lingebra.h"

#include "OCLerrors.h"
#include "OCL_device_picker.h"
#include "OCL.h"

OCLsystem cl;

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

int main(){
    int nrep = 10;
    int n    = 256;
    int n2   = n*n;

    float * A    = new float[ n2 ];
    float * B    = new float[ n2 ];
    float * C    = new float[ n2 ];
    float * Cref = new float[ n2 ];

    setArrayRandom( n2, -1.0, -1.0, A );
    setArrayRandom( n2, -1.0, -1.0, B );
    setArray      ( n2,  0.0,       Cref );
    //setArray      ( cl.buffers[2].n,  0.0,       (float*)cl.buffers[2].C_cpu );
    mmul_cpu( n2, A, B, Cref );

    cl.init();

    cl.buffers[0].n = n2;
    cl.buffers[0].typesize = sizeof(float);
    cl.buffers[0].p_cpu    = A;
    cl.buffers[0].flags    = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    cl.buffers[1].n = n2;
    cl.buffers[1].typesize = sizeof(float);
    cl.buffers[1].p_cpu    = B;
    cl.buffers[1].flags    = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;

    cl.buffers[2].n = n2;
    cl.buffers[2].typesize = sizeof(float);
    cl.buffers[2].p_cpu    = C;
    cl.buffers[2].flags    = CL_MEM_WRITE_ONLY;
    cl.buffers[2].read_on_finish = true;

    cl.initBuffers();

    cl.buildProgram( "cl/mmul/C_row.cl" );
    cl.kernels[0].kernel = clCreateKernel ( cl.program, "mmul", &cl.err );
    cl.kernels[0].dim            = 1;
    cl.kernels[0].global[0] = n2;
    cl.kernels[0].local [0] = 16;

    for (int irep = 0; irep < nrep; irep++){
        cl.kernels[0].enque( cl.commands );
    }

    cl.finish();

    printf("Hello !!! \n");

}
