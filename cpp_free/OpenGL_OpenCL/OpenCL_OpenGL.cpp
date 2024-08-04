

#include <iostream>
#include <vector>

#include <CL/opencl.hpp>
#include <GL/glew.h>
#include <SDL2/SDL.h>


const int N = 1000; // Number of bodies

// OpenGL and OpenCL variables
cl::Context clContext;
cl::CommandQueue queue;
cl::Program program;
cl::Kernel kernel_nbody;
cl::BufferGL clBuffer_pos;
cl::Buffer   clBuffer_vel;
GLuint glBuffer_pos;
GLuint glBuffer_vel;
SDL_Window* window;
SDL_GLContext glContext;

std::vector<float> positions(N*4);
std::vector<float> velocities(N*4);

void initOpenGL() {
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "Unable to initialize SDL: " << SDL_GetError() << std::endl;
        exit(1);
    }

    // Create an SDL window with OpenGL context
    window = SDL_CreateWindow("N-Body Simulation",  SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,       800, 600, SDL_WINDOW_OPENGL);
    if (!window) {   std::cerr << "Unable to create SDL window: " << SDL_GetError() << std::endl;   exit(1);  }

    // Create OpenGL context
    glContext = SDL_GL_CreateContext(window);
    if (!glContext) {    std::cerr << "Unable to create OpenGL context: " << SDL_GetError() << std::endl;  exit(1);   }

    // Initialize GLEW
    glewInit();

    // Create an OpenGL buffer
    glGenBuffers(1, &glBuffer_pos);
    glBindBuffer(GL_ARRAY_BUFFER, glBuffer_pos);
    glBufferData(GL_ARRAY_BUFFER, N * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void initOpenCL() {
    // Get all platforms (drivers)
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.empty()) {     std::cerr << "No OpenCL platforms found." << std::endl;  exit(1);   }

    // Get the first platform
    cl::Platform platform = platforms[0];

    // Get all devices
    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
    if (devices.empty()) { std::cerr << "No OpenCL devices found." << std::endl;   exit(1);  }

    // Create an OpenCL context with GL sharing
    cl_context_properties properties[] = {    CL_GL_CONTEXT_KHR, (cl_context_properties)SDL_GL_GetCurrentContext(),     CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(),    0  };
    clContext = cl::Context(devices, properties);

    queue = cl::CommandQueue(clContext, devices[0]);   // Create a command queue

    clBuffer_pos = cl::BufferGL(clContext, CL_MEM_READ_WRITE, glBuffer_pos );                      // Create a buffer from      OpenGL buffer
    clBuffer_vel = cl::Buffer(clContext, CL_MEM_READ_WRITE, N * sizeof(cl_float4), 0 );        // Create a buffer without  OpenGL buffer


    // Load and build the OpenCL program
    
    std::string kernelCode = R"(
        __kernel void nbody(__global float4* pos, __global float4* vel, float dt ){
            const int iG = get_global_id(0);
            const int n  = get_global_size(0);
            const float4 p = pos[iG];
            //const float4 v = vel[iG];
            //const float4 f = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

            if( iG==0){ printf( "CL: pos[%i](%g,%g,%g)\n", iG, p.x, p.y, p.z ); }

            // for (int j = 0; j < n; j++) {
            //     if (iG != j) {
            //         const float4 pj = pos[j];
            //         const float4 r = p - pj;
            //         const float r2 = dot(r, r);
            //         const float r3 = r2 * sqrt(r2);
            //         const float4 f_ij = r * (1.0f / r3);
            //         f += f_ij;
            //     }
            // }
            
            // v += f * dt;
            // p += v * dt;

            // if( iG==0){ printf( "CL: pos[%i](%g,%g,%g)\n", iG, p.x, p.y, p.z ); }
            
            // pos[iG] = p;
            // vel[iG] = v;
        };"
    )";
    
    program = cl::Program(clContext, kernelCode);

    //cl_int buildErr = program.build(devices);
    cl_int buildErr = program.build(devices, "-cl-std=CL1.2 -cl-mad-enable -cl-fast-relaxed-math");
    if (buildErr != CL_SUCCESS) {
        std::cerr << "Error building: " << buildErr << std::endl;
        std::string buildLog;
        program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
        std::cerr << "Build Log:" << std::endl << buildLog << std::endl;
        exit(1);
    }

    /*
    try {
        program.build(devices);
    } catch (const cl::Error &err) {
        std::cerr << "Build Error: " << err.what() << "(" << err.err() << ")" << std::endl;
        std::string buildLog;
        program.getBuildInfo(devices[0], CL_PROGRAM_BUILD_LOG, &buildLog);
        std::cerr << "Build Log:" << std::endl << buildLog << std::endl;
        exit(1);
    }
    */

    // Create the kernel
    kernel_nbody = cl::Kernel(program, "nbody");
    kernel_nbody.setArg( 0, clBuffer_pos );
    kernel_nbody.setArg( 1, clBuffer_vel );
    kernel_nbody.setArg( 2, 0.01f    );
}

void initData() {
    // Initialize positions and velocities
    for (int i=0; i<N; i++) {
        positions[i*4+0] = (float)(rand()) / RAND_MAX * 2.0f - 1.0f;
        positions[i*4+1] = (float)(rand()) / RAND_MAX * 2.0f - 1.0f;
        positions[i*4+2] = (float)(rand()) / RAND_MAX * 2.0f - 1.0f;
        positions[i*4+3] = 1.0f ;
    }
    for (int i=0; i<N; i++) { velocities[i] = 0.0f; };

    // Update OpenGL buffer with initial positions
    glBindBuffer(GL_ARRAY_BUFFER, glBuffer_pos);
    glBufferSubData(GL_ARRAY_BUFFER, 0, N * 4 * sizeof(float), positions.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Initialize OpenCL buffer with initial velocities
    queue.enqueueWriteBuffer(clBuffer_vel, CL_TRUE, 0, N * sizeof(cl_float4), velocities.data());
}

void update() {
    printf( "update() \n" );
    // Acquire the OpenGL buffer for OpenCL
    std::vector<cl::Memory> clObjects;
    clObjects.push_back(clBuffer_pos);
    //clObjects.push_back(clBuffer_vel);

    cl::Event event;
    queue.enqueueAcquireGLObjects(&clObjects);
    queue.enqueueNDRangeKernel(kernel_nbody, cl::NullRange, cl::NDRange(N), cl::NullRange);
    queue.enqueueReleaseGLObjects(&clObjects);
    queue.finish();

    cl_int status;
    event.getInfo(CL_EVENT_COMMAND_EXECUTION_STATUS, &status);
    if (status != CL_COMPLETE) {
        std::cerr << "Kernel execution failed with status: " << status << std::endl;
    }

    // Render the OpenGL buffer
    glClear(GL_COLOR_BUFFER_BIT);
    glBindBuffer(GL_ARRAY_BUFFER, glBuffer_pos);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(4, GL_FLOAT, 0, nullptr);
    glDrawArrays(GL_POINTS, 0, N);
    glDisableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    SDL_GL_SwapWindow(window);
}

int main() {
    initOpenGL();
    initOpenCL();
    initData();

    bool running = true;
    SDL_Event event;
    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
        update();
    }

    // Clean up
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
