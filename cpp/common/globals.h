#pragma once

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
        printf("  => execute action: %s", #action ); \
        {action;} \
        if(exit_on_error){ exit(0);} \
    }}

// #define _assert( pre_cond, cond, ... ){ \
//     pre_cond; \
//     if( !(cond) ){ \
//         printf("Assertion failed: %s, line: %d  function: %s file: %s \n", #cond, __LINE__, __FUNCTION__, __FILE__ ); \
//         printf("  => execute action: %s\n", #__VA_ARGS__ ); \
//         {__VA_ARGS__;} \
//         if(exit_on_error){ exit(1);} \
//     }}
#else
static bool exit_on_error = false;
//#define _assert( pre_cond, cond, action )
#define _assert( pre_cond, cond, ... )
#endif



//#endif
