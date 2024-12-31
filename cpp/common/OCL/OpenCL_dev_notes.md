## Optimize OpenCL call overhead

According to OpenCL specification Kernel argument bindings persist until:
 - The kernel object is deleted
 - The same argument index is set again
 - The program associated with kernel is recompiled

This means:
 - You don't need to set arguments again for multiple executions of the same kernel if the arguments haven't changed
 - If you switch between different kernels and come back, the previous kernel's arguments are still bound
 - Only when you want to change an argument value (like updating buffer pointers or scalar values) do you need to call clSetKernelArg again

So to optimize your code, you could:
- Set up all kernel arguments once during initialization
- Then just call clEnqueueNDRangeKernel in your main loop

Here's an example optimization pattern:
```C++
CopyInsert
// During initialization:
void init_updateJacobi_lin() {
    cl_kernel ker = kernels[iker_updateJacobi_lin];
    err |= clSetKernelArg(ker, 0, sizeof(int), &nPoint);
    err |= clSetKernelArg(ker, 1, sizeof(int), &nNeighMax);
    // ... set all static arguments ...
}

// In the main loop, just call:
void run_updateJacobi_lin_fast(int ibuff_ps_in, int ibuff_ps_out, int ibuff_bvec) {
    cl_kernel ker = kernels[iker_updateJacobi_lin];
    // Only update arguments that change:
    err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[ibuff_ps_in ].p_gpu));
    err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[ibuff_ps_out].p_gpu));
    err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_bvec  ].p_gpu));
    
    size_t local_work_size  = 32;
    size_t global_work_size = roundUp(nPoint, local_work_size);
    err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
}
```
This can significantly reduce the overhead, especially if most arguments remain constant between kernel executions.

