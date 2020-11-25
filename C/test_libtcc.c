
// see : https://stackoverflow.com/questions/9700022/tinyc-compiler-libtcc-how-to-bound-check

/*
 * Simple Test program for libtcc
 *
 * libtcc can be useful to use tcc as a "backend" for a code generator.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "libtcc.h"
 
int add(int a, int b){  return a + b; }   // this function is called by the generated code

char my_program[] =
"int fib(int n){                                  \n"
"    if (n <= 2){ return 1;                   }   \n"
"    else       { return fib(n-1) + fib(n-2); }   \n"
"}  \n"
"int foo(int n){                                  \n"
"    printf(\"Hello World!\\n\");                 \n"
"    printf(\"fib(%d) = %d\\n\", n, fib(n));      \n"
"    printf(\"add(%d, %d) = %d\\n\", n, 2 * n, add(n, 2 * n));    \n"
"    return 0;   \n"
"}  \n";

inline uint64_t getCPUticks(){
    uint32_t lo, hi;
    __asm__ __volatile__ (
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx" );
    return (uint64_t)hi << 32 | lo;
}

int main(int argc, char **argv){

    TCCState *tcc;
    int (*func)(int);
    void *mem;
    int size;

    uint64_t t1,t2;

    tcc = tcc_new();
    if (!tcc) { fprintf(stderr, "Could not create tcc state\n"); exit(1); }

    // if tcclib.h and libtcc1.a are not installed, where can we find them 
    if (argc == 2 && !memcmp(argv[1], "lib_path=",9)) tcc_set_lib_path(tcc, argv[1]+9);

    tcc_set_output_type(tcc, TCC_OUTPUT_MEMORY);               // MUST BE CALLED before any compilation

    t1 = getCPUticks();
    if (tcc_compile_string(tcc, my_program) == -1) return 1;

    // as a test, we add a symbol that the compiled program can use. You may also open a dll with tcc_add_dll() and use symbols from that 
    tcc_add_symbol(tcc, "add", add);
    size = tcc_relocate(tcc, NULL);      // get needed size of the code
    if (size == -1) return 1;
    mem = malloc(size);                  // allocate memory 
    tcc_relocate(tcc, mem);              // copy the code into it
    func = tcc_get_symbol(tcc, "foo");   // get entry symbol
    if (!func) return 1;
    tcc_delete(tcc);                     // delete the state
    t2 = getCPUticks();
    printf( "tcc compilation time: %li[ticks] %f[ms]\n", (t2-t1), (t2-t1)*1e-6 );

    t1 = getCPUticks(); 
    func(32);                            // run the code
    t2 = getCPUticks();
    printf( "run time: %li[ticks] %f[ms] \n", t2-t1, (t2-t1)*1e-6 );

    free(mem);
    return 0;
}