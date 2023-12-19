
#ifndef  OCLerrors_h
#define  OCLerrors_h

#define     OCL_DEBUG      1
//#define     OCL_DEBUG        0
static bool bOCLCheckError = true;

//----------------------------------------------------------------------------
// Purpose:  Function to output descriptions of errors for an input error code
//           and quit a program on an error with a user message
//
// RETURN:   echoes the input error code / echos user message and exits
//
// HISTORY:  Written by Tim Mattson, June 2010
//           This version automatically produced by genErrCode.py
//           script written by Tom Deakin, August 2013
//           Modified by Bruce Merry, March 2014
//           Updated by Tom Deakin, October 2014
//               Included the checkError function written by
//               James Price and Simon McIntosh-Smith
//
//----------------------------------------------------------------------------

#include <cstdio>
//#define CL_TARGET_OPENCL_VERSION 200
#include <CL/cl.h>

const char *OCL_err_code (cl_int err_in){
    switch (err_in) {
        case CL_SUCCESS:                        return (char*)"CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:               return (char*)"CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE:           return (char*)"CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:         return (char*)"CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:  return (char*)"CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:               return (char*)"CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:             return (char*)"CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:   return (char*)"CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:               return (char*)"CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:          return (char*)"CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:     return (char*)"CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:          return (char*)"CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE:                    return (char*)"CL_MAP_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:   return (char*)"CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:return (char*)"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case CL_INVALID_VALUE:                  return (char*)"CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:            return (char*)"CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:               return (char*)"CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:                 return (char*)"CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:                return (char*)"CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:       return (char*)"CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:          return (char*)"CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:               return (char*)"CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:             return (char*)"CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:return (char*)"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:             return (char*)"CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:                return (char*)"CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:                 return (char*)"CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:          return (char*)"CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:                return (char*)"CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:     return (char*)"CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:            return (char*)"CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:      return (char*)"CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:                 return (char*)"CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:              return (char*)"CL_INVALID_ARG_INDEX";
        case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:  return (char*)"CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
        case CL_INVALID_ARG_VALUE:              return (char*)"CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:               return (char*)"CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:            return (char*)"CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:         return (char*)"CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:        return (char*)"CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:         return (char*)"CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET:          return (char*)"CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:        return (char*)"CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:                  return (char*)"CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:              return (char*)"CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:              return (char*)"CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:            return (char*)"CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:              return (char*)"CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE:       return (char*)"CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_PROPERTY:               return (char*)"CL_INVALID_PROPERTY";

        default:
            return (char*)"UNKNOWN ERROR";
    }
}

void OCL_buildProgramFailure( cl_program program, cl_device_id device ){
    size_t log_size;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
    char *log = new char[log_size];
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
    printf("%s\n", log);
    delete [] log;
}

void OCL_check_error(cl_int err, const char *operation=0, char *filename=0, int line=-1, int i=-999999, const char* name=0){
    //printf( "OCL_check_error(%s) err=%i \n", operation, err );
    if (err != CL_SUCCESS){        
        //if(i==-999999){ fprintf(stderr, "Error during operation '%s'",     operation   ); }
        //else   { fprintf(stderr, "Error during operation '%s'[%i] '%s'", operation, i, name ); }
        //fprintf(stderr, "in '%s' on line %d\n", filename, line);
        //fprintf(stderr, "Error code was \"%s\" (%d)\n", OCL_err_code(err), err);
        //if(i==-999999){ fprintf(stderr, "OCL_ERROR(%i) %s during operation '%s'", err, OCL_err_code(err), operation ); }
        //else          { fprintf(stderr, "OCL_ERROR(%i) %s during operation '%s'", err, OCL_err_code(err), operation, i, name );  }
        fprintf(stderr, "OCL_ERROR(%i) %s", err, OCL_err_code(err));
        if(operation)fprintf(stderr, " during %s", operation );
        if(filename )fprintf(stderr, "%s #line=%i", filename, line );
        if(name     )fprintf(stderr, "%s[%i]", name, i );
        fprintf(stderr, "\n" );
        exit(0);
    }
}

void OCL_error_warn(cl_int err, const char *operation, int i=-999999 ){
    if (err != CL_SUCCESS){        
        if(i==-999999){ fprintf(stderr, "WARRNING: Error in '%s'     \"%s\" (%d)\n", operation   , OCL_err_code(err), err ); }
        else          { fprintf(stderr, "WARRNING: Error in '%s'[%i] \"%s\" (%d)\n", operation, i, OCL_err_code(err), err ); }
    }
}

void printBinary( int err ){
    for(int i=0;i<32;i++){
        err>>=1;
        printf("%1i", err&1 );
    }
    printf("\n" );
}

#if OCL_DEBUG
#define OCLerr(E)                                       OCL_check_error(E,"",__FILE__,__LINE__);
#define OCL_checkError(E, S)          if(bOCLCheckError)OCL_check_error(E,S,__FILE__,__LINE__);
#define OCL_checkError_(E, S,I)       if(bOCLCheckError)OCL_check_error(E,S,__FILE__,__LINE__,I);
#define OCL_checkError__(E, S,I,name) if(bOCLCheckError)OCL_check_error(E,S,__FILE__,__LINE__,I,name);
#define OCL_error_warn_(E,S,I)        if(bOCLCheckError)OCL_error_warn(E,S,I);

#else
#define OCLerr(E)               {}
#define OCL_checkError(E, S)    {}
#define OCL_checkError_(E, S,I) {}
#define OCL_checkError__(E, S,I,name) {}
#endif

#endif
