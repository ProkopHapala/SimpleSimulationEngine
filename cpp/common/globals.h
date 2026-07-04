#pragma once

/// @file globals.h
/// @brief Global variables and debug utilities shared across all translation units.
///
/// This file is included everywhere — it holds cross-cutting state that doesn't belong to any
/// single module: verbosity/idebug levels (controlled at runtime to gate debug output), tmpstr
/// (shared formatting buffer to avoid per-call stack allocation), tick2second (hardware timer
/// calibration), and DEBUG/DBG print macros.
///
/// The most important utility here is _assert(pre_cond, cond, action) — a three-argument
/// assertion macro designed for scientific computing where error context is essential:
/// - pre_cond runs before the check (e.g. `int iv=findVert(pos)` — does the lookup)
/// - cond is the actual assertion (e.g. `iv<0` — vertex should NOT exist)
/// - action runs on failure (e.g. diagnostic printf with both positions and distance)
///
/// In DEBUGBUILD, failure prints location info, executes action, and exits if exit_on_error.
/// In release, the entire macro (including pre_cond) is a no-op — so pre_cond must NOT have
/// side effects that the program depends on. This is by design: debug checks should not
/// affect release behavior.

//#ifndef  globals_h

//extern 
static int  verbosity = 1;
static int  idebug    = 0;

constexpr static const int ntmpstr=1024;
static char tmpstr[ntmpstr];

static double tick2second=1e-9;

#define DEBUG   printf( "DEBUG #l %i %s \n",    __LINE__, __FUNCTION__ );
#define DEBUGF  printf( "DEBUG #l %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ );
#define DBG(format,args...) { printf("DEBUG "); printf(format,## args); }


// depending on debug / optimization leval 
#ifdef DEBUGBUILD
static bool exit_on_error = true;

#define _assert( pre_cond, cond, action ){ \
    pre_cond; \
    if( !(cond) ){ \
        printf("Assertion failed: %s, line: %d  function: %s file: %s \n", #cond, __LINE__, __FUNCTION__, __FILE__ ); \
        printf("  => execute action:\n %s\n", #action ); \
        {action;} \
        if(exit_on_error){ exit(0);} \
    }}

// #define _assert( pre_cond, cond, ... ){ \
//     pre_cond; \
//     if( !(cond) ){ \
//         printf("Assertion failed: %s, line: %d  function: %s file: %s \n", #cond, __LINE__, __FUNCTION__, __FILE__ ); \
//         printf("  => execute action\n %s\n", #__VA_ARGS__ ); \
//         {__VA_ARGS__;} \
//         if(exit_on_error){ exit(1);} \
//     }}
#else
static bool exit_on_error = false;
//#define _assert( pre_cond, cond, action )
#define _assert( pre_cond, cond, ... )
#endif



//#endif
