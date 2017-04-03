#ifndef  convolve2D_h
#define  convolve2D_h

//============== Globals

constexpr int nx   = 256+2;
constexpr int ny   = 256+2;
constexpr int ntot = nx*ny;

float buff [nx*ny];
float buff_[nx*ny];

void genXOR2D(int nx, int ny, float * data){
    float scale_val = 1.0/256.0;
    for(int iy=0;iy<ny;iy++){ for(int ix=0;ix<nx;ix++){
        buff[iy*nx+ix] = ((ix^iy)&0xFF)*scale_val;
        //data[iy*nx+ix] = 0.5;
    } }
}

void blur(int nx, int ny, float * I, float * O ){
    float renorm = 1.0/9.0;
    for(int iy=1;iy<ny-1;iy++){ for(int ix=1;ix<nx-1;ix++){
        int i   = iy*nx+ix;
        O[i] =( I[i-nx-1] + I[i-nx] + I[i-nx+1] +
                I[i   -1] + I[i   ] + I[i   +1] +
                I[i+nx-1] + I[i+nx] + I[i+nx+1] ) * renorm;
    } }
}

void blur_Gauss(int nx, int ny, float * I, float * O ){
    for(int iy=1;iy<ny-1;iy++){ for(int ix=1;ix<nx-1;ix++){
        int i   = iy*nx+ix;
        //O[i] =( 0.0625*I[i-nx-1] + 0.125*I[i-nx] + 0.0625*I[i-nx+1] +
        //        0.125 *I[i   -1] + 0.25 *I[i   ] + 0.125 *I[i   +1] +
        //        0.0625*I[i+nx-1] + 0.125*I[i+nx] + 0.0625*I[i+nx+1] );

        O[i] = 0.25  *  I[i]
         + 0.0625*( I[i-nx-1] + I[i+nx-1] + I[i+nx+1] + I[i-nx+1] )
         + 0.125 *( I[i-nx]   + I[i-1]    + I[i+1]    + I[i+nx]   );
    } }
}


//============== Functions

#endif

