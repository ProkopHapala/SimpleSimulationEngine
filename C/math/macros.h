#ifndef  macros_h
#define  macros_h

#define _inline static inline

// see https://stackoverflow.com/questions/1253934/c-pre-processor-defining-for-generated-function-names
//#define MAKE_FN_NAME(x) void  Callback_ ## x (void)
//#define FUNCTION_NAME(signal) MAKE_FN_NAME(signal)
#define METHOD_CONCAT(TYPE,rest)  TYPE ## _ ## rest
#define METHOD(TYPE,rest)         METHOD_CONCAT(TYPE,rest)
//#define METHOD_NAME_CONCAT(TYPE,name)  TYPE ## _ ## name
//#define METHOD_NAME(TYPE,name)         METHOD_NAME_CONCAT(TYPE,name)
//#define METHOD(TYPE,rest) TYPE ## _ ## rest 

//#MANGLE_FUNC (pre,rest) pre VEC rest

#endif