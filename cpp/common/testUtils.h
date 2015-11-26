
/*
void printVec( const Vec3f& vec ){
	printf( " %f %f %f ", vec.x, vec.y, vec.z );
}

void printMat( const Mat3f& mat  ){
	printf( " %f %f %f \n", mat.xx, mat.xy, mat.xz );
	printf( " %f %f %f \n", mat.yx, mat.yy, mat.yz );
	printf( " %f %f %f \n", mat.zx, mat.zy, mat.zz );
}
*/

void printVec( const Vec3d& vec ){
	printf( " %f %f %f ", vec.x, vec.y, vec.z );
}

void printQuat( const Quat4d& q ){
	printf( " %f %f %f %f ", q.x, q.y, q.z, q.w );
}

void printMat( const Mat3d& mat  ){
	printf( " %f %f %f \n", mat.ax, mat.ay, mat.az );
	printf( " %f %f %f \n", mat.bx, mat.by, mat.bz );
	printf( " %f %f %f \n", mat.cx, mat.cy, mat.cz );
}

// CPU ticks timer
// http://stackoverflow.com/questions/6432669/variance-in-rdtsc-overhead

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






