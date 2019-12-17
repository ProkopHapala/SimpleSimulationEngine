#ifndef  PixMapFractal_h
#define  PixMapFractal_h

int patterns_1[8] = {
1,0,
0,0,

0,1,
1,1
};

inline int binaryPixMapFrac( int nbit, int nbitend, int ix, int iy, int * patterns ){
    int ipatern = 0;
    for(int ibit=nbit; ibit>nbitend; ibit>>=1 ){
        int dx = (ix&ibit)>0;
        int dy = (iy&ibit)>0;
        int ipix = (ipatern<<2) + (dy<<1) + dx;
        //printf( " %i : %i (%i,%i) -> %i \n", ibit, dx,dy, ipatern, ipix );
        ipatern = patterns[ ipix ];
        //ipatern += dx*dy;
    }
    return ipatern;
}

#endif
