import os
import re
import numpy as np
import pyopencl as cl

def print_devices(platforms=None ):
    if platforms is None:
        platforms = cl.get_platforms()
    for i, platform in enumerate(platforms):
        print(f"Platform {i}: {platform.name}")
        devices = platform.get_devices()
        for j, device in enumerate(devices):
            print(f"  Device {j}: {device.name}")

def select_device(platforms=None, preferred_vendor='nvidia', bPrint=False, device_index=0):
    if platforms is None:
        platforms = cl.get_platforms()
    if bPrint:
        print_devices(platforms)
    # Try to find preferred vendor device
    preferred_devices = []
    for platform in platforms:
        for device in platform.get_devices():
            if preferred_vendor.lower() in device.name.lower():
                preferred_devices.append((platform, device))
    if preferred_devices:
        platform, device = preferred_devices[0]
        ctx = cl.Context([device])
        if bPrint:
            print(f"Selected {preferred_vendor} device: {device.name}")
    else:
        # Fall back to default behavior
        if bPrint:
            print(f"Selected default device {device_index}")
        ctx = cl.create_some_context(answers=[device_index])
    return ctx

class OCLSystem:
    """
    A system for managing OpenCL context, programs, kernels, and buffers.
    """
    def __init__(self, nloc=32, device_index=0, preferred_vendor='nvidia'):
        self.nloc = nloc
        self.ctx = select_device(preferred_vendor=preferred_vendor, bPrint=True, device_index=device_index)
        self.queue = cl.CommandQueue(self.ctx)
        self.programs = {}
        self.buffers = {}
        self.kernels = {} # Global cache for kernel objects by name
        self.kernel_params = {}

    def load_program(self, name, kernel_filepath):
        """
        Load, compile, and cache all kernels from an OpenCL program file.
        Args:
            name (str): A unique name for this program.
            kernel_filepath (str): Absolute path to the kernel file.
        """
        print(f"OCLSystem::load_program() Loading kernel program '{name}' from: {kernel_filepath}")
        if not os.path.exists(kernel_filepath):
            raise FileNotFoundError(f"OCLSystem::load_program() ERROR: Kernel file not found at: {kernel_filepath}")
        with open(kernel_filepath, 'r') as f:
            kernel_source = f.read()
            try:
                prg = cl.Program(self.ctx, kernel_source).build()
                self.programs[name] = prg
                # NEW: Cache all kernels from this program globally by their name
                for kernel_func in prg.all_kernels():
                    kernel_name = kernel_func.function_name
                    if kernel_name in self.kernels:
                        print(f"Warning: Kernel '{kernel_name}' is being redefined. The last loaded version will be used.")
                    self.kernels[kernel_name] = kernel_func
                print(f"OCLSystem::load_program() Successfully loaded and cached kernels from program '{name}'.")
                return True
            except cl.LogicError as e:
                print(f"OCLSystem::load_program() OpenCL compilation error for '{name}': {e}")
                raise

    def create_buffer(self, name, size, flags=cl.mem_flags.READ_WRITE, hostbuf=None):
        """
        Create an OpenCL buffer and store it.
        """
        if size <= 0:
            self.buffers[name] = None
            return None
        
        if name in self.buffers and self.buffers[name] is not None:
            if self.buffers[name].size != size:
                print(f"OCLSystem::create_buffer() Re-allocating buffer '{name}' with new size {size} bytes (old size {self.buffers[name].size} bytes)")
                self.buffers[name].release()
                self.buffers[name] = cl.Buffer(self.ctx, flags, size=size, hostbuf=hostbuf)
        else:
            print(f"OCLSystem::create_buffer() Allocating buffer '{name}' with size {size} bytes")
            self.buffers[name] = cl.Buffer(self.ctx, flags, size=size, hostbuf=hostbuf)
        return self.buffers[name]

    def toGPU(self, buffer_name, host_data, byte_offset=0):
        """
        Upload data to a GPU buffer.
        """
        buffer_obj = self.buffers.get(buffer_name)
        if buffer_obj is None:
            raise ValueError(f"Buffer '{buffer_name}' does not exist or is zero-sized.")
        cl.enqueue_copy(self.queue, buffer_obj, host_data, device_offset=byte_offset).wait()

    def fromGPU(self, buffer_name, host_data, byte_offset=0):
        """
        Download data from a GPU buffer.
        """
        buffer_obj = self.buffers.get(buffer_name)
        if buffer_obj is None:
            raise ValueError(f"Buffer '{buffer_name}' does not exist or is zero-sized.")
        cl.enqueue_copy(self.queue, host_data, buffer_obj, device_offset=byte_offset).wait()

    def get_buffer(self, name):
        """
        Get a buffer by name.
        """
        return self.buffers.get(name)

    def get_kernel(self, name):
        """
        Get a compiled kernel object by its function name from the global cache.
        """
        kernel_obj = self.kernels.get(name)
        if kernel_obj is None:
            raise ValueError(f"Kernel function '{name}' not found in any loaded program.")
        return kernel_obj

    def execute_kernel(self, kernel_obj, global_size, local_size, *args):
        """
        Execute a pre-fetched OpenCL kernel object.
        """
        kernel_obj(self.queue, global_size, local_size, *args).wait()

    def set_kernel_param(self, name, value):
        """
        Set a scalar parameter for kernels.
        """
        self.kernel_params[name] = value

    def get_kernel_param(self, name):
        """
        Get a scalar parameter by name.
        """
        return self.kernel_params.get(name)

    def clear(self):
        """
        Clear all allocated resources.
        """
        print("OCLSystem::clear() Releasing all OpenCL resources.")
        self.buffers.clear()
        self.kernel_params.clear()
        self.kernels.clear()
        self.programs.clear()