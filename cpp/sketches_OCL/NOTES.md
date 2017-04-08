
# Notes

- atomic operations on GPU
  http://simpleopencl.blogspot.cz/2013/04/performance-of-atomics-atomics-in.html
- dventures in OpenCL Part 3: Constant Memory Structs
http://enja.org/2011/03/30/adventures-in-opencl-part-3-constant-memory-structs/
- bank conflicts
http://stackoverflow.com/questions/3841877/what-is-a-bank-conflict-doing-cuda-opencl-programming
- async_work_group_copy
    - http://stackoverflow.com/questions/43217449/shouldnt-be-3x3-convolution-much-faster-on-gpu-opencl
    - https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/async_work_group_copy.html
    - http://stackoverflow.com/questions/15545841/how-to-use-async-work-group-copy-in-opencl
    - https://streamcomputing.eu/blog/2014-06-19/using-async_work_group_copy-on-2d-data/

# Examples & Tutorials

 - http://nullprogram.com/tags/gpgpu/
    2014-06-29 A GPU Approach to Particle Physics
    2014-06-22 A GPU Approach to Path Finding
    2014-06-10 A GPU Approach to Conway's Game of Life
    2014-06-01 A GPU Approach to Voronoi Diagrams 
 -   Brief Intro to GPU PIC with CUDA particle in cell GPU https://www.particleincell.com/2016/cuda-pic/\
 - FFT OpenCL
http://www.bealto.com/gpu-fft.html
 
# Libraries

- https://github.com/ComputationalRadiationPhysics/picongpu
  PIConGPU is a fully relativistic, many GPGPU, 3D3V particle-in-cell (PIC) code. The Particle-in-Cell algorithm is a central tool in plasma physics. It describes the dynamics of a plasma by computing the motion of electrons and ions in the plasma based on Maxwell's equations.
  - opencl particle in cell on nvidia examples
https://github.com/sschaetz/nvidia-opencl-examples/blob/master/OpenCL/src/oclParticles/Particles.cl
  - http://scicomp.stackexchange.com/questions/24113/mixing-some-particles-together-game-physics-for-engineers 
  - GPU Molecular dynamics
	https://github.com/pandegroup/openmm
	https://github.com/Mantevo/miniMD
	https://github.com/NVIDIA/CoMD-CUDA
	https://github.com/AccelerateHS/accelerate

# Resources & References

 - FREE HTML version
The OpenCL Programming Book. Copyright © 2010 Fixstars Corporation. All rights reserved.
https://www.fixstars.com/en/opencl/book/OpenCLProgrammingBook/contents/
free-programming-books
https://github.com/vhf/free-programming-books/blob/master/free-programming-books.md
- OpenCL Resources
https://github.com/KhronosGroup/Khronosdotorg/blob/master/api/opencl/resources.md
- non uniform histogram implementation (not OpenCL)
http://stackoverflow.com/questions/15707064/efficient-histogram-implementation-using-a-hash-function
- Lattice boltzman OpenCL
http://scicomp.stackexchange.com/questions/19740/gpu-enabled-lattice-boltzmann-solvers
https://github.com/sailfish-team/sailfish

## Tasks:
  - 3x3 convolution
    - https://software.intel.com/en-us/blogs/2014/07/15/an-investigation-of-fast-real-time-gpu-based-image-blur-algorithms
     - http://gpgpu2.blogspot.cz/
     - http://downloads.ti.com/mctools/esd/docs/opencl/optimization/examples.html
     - http://stackoverflow.com/questions/43217449/shouldnt-be-3x3-convolution-much-faster-on-gpu-opencl/43217833?noredirect=1#comment73553552_43217833


## Unsorted & Temp






