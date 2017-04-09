
// https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/read_imagef2d.html
// https://github.com/codehathi/CUDA-OpenCL-Hardware-Interpolation
// http://www.cmsoft.com.br/opencl-tutorial/opencl-image2d-variables/
// http://amdahlsoftware.com/matrix-multiplication-using-the-opencl-image-type/

// /home/prokop/git_SW/_OpenCL/CUDA-OpenCL-Hardware-Interpolation/OpenCL/src
// /home/prokop/Dropbox/MyDevSW/OpenCL/OpenCL_in_C/MMUL_images

// copy buffer to image   https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clEnqueueCopyBufferToImage.html
// usefull especially if direct read-write 3D image is not supported (in OpenCL 1.0)   https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/cl_khr_3d_image_writes.html

__constant sampler_t s0 = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

__kernel void getValInPoints(
    __read_only image2d_t  hf,
    __global  float2*      points,
    __global  float *      vals
){
    const float2 coord     = points[get_global_id(0)];
    vals[get_global_id(0)] = read_imagef(hf, s0, coord).w;
}

__kernel void getF2InPoints(
    __read_only image2d_t  hf,
    __global  float2*      points,
    __global  float2*      Dvals
){
    const float2 coord      = points[get_global_id(0)];
    Dvals[get_global_id(0)] = read_imagef(hf, s0, coord).xy;
}

__kernel void relaxPoints(
    __read_only image2d_t  hf,
    __global  float2*      points,
    __global  float2*      Dvals,
    int niters, float dt, float damp
){
    //if(get_global_id(0)==0) printf("relaxPoint %i %f %f \n", niters, damp, dt );
    float2 p  = points[get_global_id(0)];
    float2 v  = (float2)(0.0f,0.0f);
    float2 f; 
    for(int i=0; i<niters; i++){
        f  = read_imagef(hf, s0, p).xy;
        v *= damp;
        v -= f*dt;
        p += v*dt;
    }
    points[get_global_id(0)] = p;
    //points[get_global_id(0)] = 0;
    //Dvals[get_global_id(0)] = f;
}

