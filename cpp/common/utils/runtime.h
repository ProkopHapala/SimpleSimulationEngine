
#ifndef runtime_h
#define runtime_h

#include "IO_utils.h"

// see http://www.yolinux.com/TUTORIALS/LibraryArchives-StaticAndDynamic.html

void recompileLib( const char * dbin, const char * dsrc, const char * fname, const char * fflags ){
    //system ( "echo \"recompilation START\"; $(cd plugin; make); echo \"recompilation DONE\";" );
    printf("Compiling '%s' | dbin '%s'  dsrc '%s' fflags '%s' \n", fname, dbin, dsrc, fflags );
    char buf[8096];
    sprintf( buf,
//     "cd plugins/;"
     "cd %s;"
     "rm lib%s.so;"
     "rm %s.o;"
//     "gcc -Wall -fPIC -c %s.c;"
     "g++ -std=c++11 %s -fPIC -c %s%s.cpp;"
     "g++ -shared -Wl,-soname,lib%s.so -o lib%s.so %s.o;"
     , dbin, fname, fname,   fflags, dsrc, fname,     fname, fname, fname
    );
    printf("Compilation string >>>\n %s \n<<<\n", buf );
    printf("==== BEGIN COMPILATION ==== \n");
    system( buf );
    printf("==== END   COMPILATION ==== \n");
}

#endif

