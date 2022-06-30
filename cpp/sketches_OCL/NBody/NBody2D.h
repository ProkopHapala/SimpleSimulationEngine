#ifndef  NBody2D_h
#define  NBody2D_h

#include "arrayAlgs.h"

#define R2SAFE  1.0e-2
#define F2MAX   10.0

//  ============== Globals

constexpr float R_MAX = 1.8;
constexpr float R2MAX = R_MAX*R_MAX;

//int n = 1024;
//constexpr int n = 64;
//constexpr int n = 128;
//constexpr int n = 256;
//constexpr int n = 512;
//constexpr int n = 1024;
//constexpr int n = 2048;
//constexpr int n = 4096;
//constexpr int   n = 8192;
constexpr int n = 16384;
//constexpr int n = 32768;
float time_step = 0.05;
//constexpr int n = 1024*2;
float damp      = 0.5;
Vec2f pos   [n];
Vec2f pos_  [n];
Vec2f vel   [n];
Vec2f force [n];
Vec2f force_[n];

// cell acceleration buffer
constexpr int  nx    = 64;
constexpr int  ny    = 64;

//constexpr int  nx    = 8;
//constexpr int  ny    = 8;
constexpr int ncell = nx*ny;
constexpr float cell_size = 4.0;
constexpr float xspan = cell_size*(ny-2)*0.49;
constexpr float yspan = cell_size*(ny-2)*0.49;

int           cell2pos  [ncell+1];
int           cellNs    [ncell];

//int           cellNcount[ncell];
int           atom2cell[n];

float grid_orig [2] = {-cell_size*nx*0.5,-cell_size*ny*0.5};

//constexpr int nbins           = 4;

constexpr int nbins               = 8;
const     int binTrashs[nbins]    = {4,8,16,32,64,128,256};
const     int cellBinTrash[nbins] = {32,64,128,256,512,1024,2048,4096};
int           istrata[nbins];
std::vector<int>   cellBins[nbins];
std::vector<Vec2i> interBins[nbins];

int   ninters=0;
int   interNs    [ncell*9];
Vec2i interIJs   [ncell*9];
int   interPermut[ncell*9];

int   sortedCells[ncell];
int   sortedIters[ncell*4*9];

#include <map>
std::multimap<int,Vec2i> interMap;

const int neighCells[9] = {
    0,    -1,    +1,
    -nx-1,-nx,-nx+1,
    +nx-1,+nx,+nx+1
};

int       nFilledCells = 0;
Vec2i     cellBounds     [ncell];
//int       DEBUG_int_buff [ncell];
//int       DEBUG_int_buff_[ncell];
//int DEBUG_counter = 0;

//  ============== Functions

inline int getBinIndex(int i){
    //if i<binTrashs[1]{ if(i<binTrashs[0]){ return 0; }else{ return 1; }
    //}else{             if(i<binTrashs[2]){ return 2; }else{ return 3; }   };
    if(i<32){
        if (i<8)  { if(i<4 ){ return 0; }else{ return 1; }
        }else     { if(i<16){ return 2; }else{ return 3; }   };
    }else{
        if (i<128){ if(i<64 ){ return 4; }else{ return 5; }
        }else     { if(i<256){ return 6; }else{ return 7; }   };
    }
}


inline int getCellBinIndex(int i){
    //if i<binTrashs[1]{ if(i<binTrashs[0]){ return 0; }else{ return 1; }
    //}else{             if(i<binTrashs[2]){ return 2; }else{ return 3; }   };
    if(i<cellBinTrash[3]){
        if (i<cellBinTrash[1]){ if(i<cellBinTrash[0] ){ return 0; }else{ return 1; }
        }else                 { if(i<cellBinTrash[2] ){ return 2; }else{ return 3; }   };
    }else{
        if (i<cellBinTrash[5]){ if(i<cellBinTrash[4] ){ return 4; }else{ return 5; }
        }else                 { if(i<cellBinTrash[6] ){ return 6; }else{ return 7; }   };
    }
}


void sizeSort( ){
    printf(" --- sizeSort \n");
    // http://stackoverflow.com/questions/53161/find-the-highest-order-bit-in-c
    for(int i=0; i<nbins; i++ ){ cellBins[i].clear(); }
    int nsum = 0;
    for(int i=0; i<ncell; i++ ){
        int ni = cellNs[i];
        if(ni==0) continue;
        nsum+=ni;
        //int firstBit = fls(ibin);
        int ibin = getBinIndex(ni);
        cellBins[ibin].push_back(i);
        //printf("%i %i %i\n", i, ni, ibin);
    }
    printf("nsum %i\n", nsum);
    for(int i=0; i<nbins; i++ ){ printf( "nbin %i size %i (%i)\n", i, cellBins[i].size(), cellBins[i].size()*binTrashs[i]  ); }
    //exit(0);
}

void sortCells( ){
    //printf(" --- sortInterCell \n");
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    //for(int i=0; i<nbins; i++ ){ cellBins[i].clear(); }
    //interMap.clear();
    for(int ibin=0; ibin<nbins; ibin++){ cellBins[ibin].clear(); }
    //FILE * logFile = fopen("sortCells.log","w");
    for(int iy=1; iy<ny_; iy++){
        for(int ix=1; ix<nx_; ix++){
            int i  = nx*iy+ix;
            int ni = cellNs[i];
            if( ni==0 ) continue;
            int nj = ni;
            nj += cellNs[i-nx-1]; nj += cellNs[i-nx  ]; nj += cellNs[i-nx+1];
            nj += cellNs[i-1   ];                       nj += cellNs[i+1   ];
            nj += cellNs[i+nx-1]; nj += cellNs[i+nx  ]; nj += cellNs[i+nx+1];
            int ibin = getCellBinIndex(nj*ni);
            cellBins[ibin].push_back(i);
            //fprintf( logFile, "cell(%i,%i) %i (%i,%i) %i %i \n", iy, ix, i, ni, nj, ni*nj, ibin );
        }
    }
    nFilledCells = 0;
    for(int ibin=0; ibin<nbins; ibin++){
        std::vector<int>& bin = cellBins[ibin];
        for(int i=0; i<bin.size(); i++){
            int icell = bin[i];
            sortedCells[nFilledCells]=icell;
            //cellBounds[nFilledCells].set( cell2pos[icell], cellNs[icell] );
            int nj = cellNs[icell] + cellNs[icell-1] + cellNs[icell+1] + cellNs[icell-nx-1] + cellNs[icell-nx] + cellNs[icell-nx+1] + cellNs[icell+nx-1] + cellNs[icell+nx] + cellNs[icell+nx+1];
            //fprintf( logFile, "cellBounds %i %i  %i (%i,%i) %i %i \n", nFilledCells, icell, cell2pos[icell], cellNs[icell], nj, nj*cellNs[icell], ibin );
            nFilledCells++;
        }
        istrata[ibin] = nFilledCells;
    }
    //fclose(logFile);
    //printf("nsum %i\n", nsum);
    //for(int i=0; i<nbins; i++ ){ printf( "nbin %i size %i \n", i, cellBins[i].size()  ); }
    //exit(0);
}












inline void processIterCell(int i, int ni, int j ){
    int nj = cellNs[j];
    if (nj==0) return;
    int ibin = getBinIndex(ni*nj);
    interBins[ibin].push_back({i,j});
    //interMap.insert( std::pair<int,Vec2i>(ni*nj,{i,j}) );   // too costly - 16 Mticks for 4096 particles and 32*32 cells
}

void sortIntersBins( ){
    //printf(" --- sortInterCell \n");
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    //for(int i=0; i<nbins; i++ ){ cellBins[i].clear(); }
    //interMap.clear();
    for(int ibin=0; ibin<nbins; ibin++){ interBins[ibin].clear(); }
    for(int iy=1; iy<ny_; iy++){
        for(int ix=1; ix<nx_; ix++){
            int i  = nx*iy+ix;
            int ni = cellNs[i];
            if( ni==0 ) continue;
            processIterCell( i, ni, i-nx-1 );
            processIterCell( i, ni, i-nx   );
            processIterCell( i, ni, i-nx+1 );
            processIterCell( i, ni, i   -1 );
            processIterCell( i, ni, i      );
            processIterCell( i, ni, i   +1 );
            processIterCell( i, ni, i+nx-1 );
            processIterCell( i, ni, i+nx   );
            processIterCell( i, ni, i+nx+1 );
        }
    }
    //FILE * logFile = fopen("sortIntersBins.log","w");
    ninters = 0;
    for(int ibin=0; ibin<nbins; ibin++){
        std::vector<Vec2i>& bin = interBins[ibin];
        for(int i=0; i<bin.size(); i++){
            Vec2i c = bin[i];
            int i4 = ninters<<2;
            sortedIters[i4  ] = cell2pos[c.a];
            sortedIters[i4+1] = cellNs  [c.a];
            sortedIters[i4+2] = cell2pos[c.b];
            sortedIters[i4+3] = cellNs  [c.b];
            //fprintf( logFile, "inter %i (%i,%i) (%i,%i) \n", ninters, sortedIters[i4  ], sortedIters[i4+1], sortedIters[i4+2], sortedIters[i4+3] );
            ninters++;
        }
    }
    //fclose(logFile);
    //printf("nsum %i\n", nsum);
    //for(int i=0; i<nbins; i++ ){ printf( "nbin %i size %i \n", i, cellBins[i].size()  ); }
    //exit(0);
}

inline void processIterCell_2( int i, int ni, int j ){
    int nj = cellNs[i];
    if (nj==0) return;
    interNs [ninters] = ni*nj;
    interIJs[ninters] = {i,j};
    ninters++;
}

void sortIntersQuick( ){
    ninters = 0;
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    for(int iy=1; iy<ny_; iy++){
        for(int ix=1; ix<nx_; ix++){
            int i  = nx*iy+ix;
            int ni = cellNs[i];
            if( ni==0 ) continue;
            processIterCell_2( i, ni, i-nx-1 );
            processIterCell_2( i, ni, i-nx   );
            processIterCell_2( i, ni, i-nx+1 );
            processIterCell_2( i, ni, i   -1 );
            processIterCell_2( i, ni, i      );
            processIterCell_2( i, ni, i   +1 );
            processIterCell_2( i, ni, i+nx-1 );
            processIterCell_2( i, ni, i+nx   );
            processIterCell_2( i, ni, i+nx+1 );
        }
    }
    for(int i=0; i<ninters; i++ ){
        interPermut[i]=i;
        //printf( "inter-cell %i (%i,%i) %i \n", i, interIJs[i].a, interIJs[i].b, interNs[i] );
    }
    //quickSort<int>( interNs, interPermut, 0, ninters ); // too costly - 16 Mticks for 4096 particles and 32*32 cells
    //for(int i=0; i<ninters; i++ ){ int ip=interPermut[i]; printf( "inter-cell %i -> %i (%i,%i) %i \n", i, ip, interIJs[ip].a, interIJs[ip].b, interNs[ip] ); }
    //exit(0);
}



/*
inline void acum_force( const Vec2f& p1, const Vec2f& p2, Vec2f& f ){
    Vec2f dp; dp.set_sub( p2, p1 );
    float ir2 = 1/( dp.norm2() + R2SAFE );
    //float ir  = sqrt(ir2);
    float ir6 = ir2*ir2*ir2;
    float fr  = ir6 - ir6*ir6;
    f.add_mul( dp, fr );
}
*/

inline void acum_force( const Vec2f& p1, const Vec2f& p2, Vec2f& f ){
    Vec2f dp; dp.set_sub( p2, p1 );
    float r2 = dp.norm2();
    if (r2 > R2MAX) return;
    //Draw2D::drawLine(p1,p2);

    float ir2 = 1/( r2 + R2SAFE );
    //float fr  = (0.2-ir2)*(R2MAX-r2);
    float fr  = (0.7-ir2)*(R2MAX-r2);
    f.add_mul( dp, fr );
}

void initParticles( double step, double drnd ){
    int ny    = (int)(round(sqrt(n)));
    int nx    = n / ny;
    //int nrest = n - ny*nx;
    ny +=1;
    int i = 0;
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            if(i<n){
                pos[i].set(
                    (ix-0.5*nx)*step + randf(-drnd,drnd),
                    (iy-0.5*ny)*step + randf(-drnd,drnd)
                );
                vel  [i].set(0.0);
                force[i].set(0.0);
            }
            i++;
        }
    }

    //cellBins[] = std::vector<int>[nbins];
}

void atomsToCells( ){
    double inv_size = 1/cell_size;
    for(int i=0; i<ncell; i++){ cellNs[i] = 0; }
    for(int i=0; i<n; i++){
        Vec2f& p  = pos[i];
        int ix    = (p.x - grid_orig[0])*inv_size;
        int iy    = (p.y - grid_orig[1])*inv_size;
        int icell = nx*iy + ix;
        cellNs[icell]++;
        atom2cell[i] = icell;
        //printf( "%i (%3.3f,%3.3f) (%i,%i) %i \n", i, p.x,p.y, ix, iy, icell );
    }
    int nsum = 0;
    for(int i=0; i<ncell; i++){
        cell2pos[i] = nsum;
        nsum       += cellNs[i];
        cellNs[i]   = 0;
    }
    cell2pos[ncell] = n;
    //for(int i=0; i<ncell; i++){ printf( "cell2pos[%i] = %i \n", i, cell2pos[i] ); }
    for(int i=0; i<n; i++){
        int icell     = atom2cell[i];
        int ni        = cellNs  [icell];
        int i_        = cell2pos[icell] + ni;
        cellNs[icell] = ni+1;
        pos_[i_]      = pos[i];
    }
    for(int i=0; i<n; i++){ pos[i]=pos_[i]; }
}

void prepareCellBounds( ){
    for(int i=0; i<ncell; i++){ cellBounds[i].set( cell2pos[i], cellNs[i] ); }
};

void add_ineraction_forces( int n, Vec2f* pos, Vec2f* force ){
    for(int i=0; i<n; i++){
        Vec2f f; f.set(0.0);
        Vec2f& pi = pos[i];
        for(int j=0; j<n; j++){ acum_force( pi, pos[j], f ); }
        force[i] = f;
    }
}

void add_forces_Inters( int n, Vec2f* pos, Vec2f* force ){
    Vec2f pjs[256];
    //FILE * logFile = fopen("add_forces_Inters.log","w");
    for(int inter=0; inter<ninters; inter++){
        int i4 = inter<<2;
        int i0 = sortedIters[i4  ];
        int ni = sortedIters[i4+1];
        int j0 = sortedIters[i4+2];
        int nj = sortedIters[i4+3];
        //fprintf( logFile, "inter %i (%i,%i) (%i,%i) %i \n", inter, i0, j0, ni, nj, ni*nj );
        for(int j=0;j<nj;j++){ pjs[j] = pos[j0+j]; }
        for(int i=0;i<ni;i++){
            Vec2f pi = pos[i0+i];
            Vec2f f; f.set(0.0);
            for(int j=0;j<nj;j++){ acum_force( pi, pjs[j], f ); }
            //for(int j=0;j<nj;j++){ acum_force( pi, pos[j0+j], f ); }
            force[i0+i].add( f );
        }
    }
    //fclose(logFile);
}

void interact_cell( int iStart, int iEnd, int jStart, int jEnd, Vec2f* pos, Vec2f* force ){
    //DEBUG_counter += (jEnd-jStart);
    for(int i=iStart; i<iEnd; i++){
        Vec2f& pi = pos[i];
        Vec2f  f; f.set(0.0);
        for(int j=jStart; j<jEnd; j++){ acum_force( pi, pos[j], f ); }
        force[i].add(f);
    }
}

void interact_cell( int iStart, int iEnd, Vec2f* pos, Vec2f* force ){
    for(int i=iStart; i<iEnd; i++){
        Vec2f& pi = pos[i];
        Vec2f  f; f.set(0.0);
        for(int j=iStart; j<iEnd; j++){
            if(i!=j) acum_force( pi, pos[j], f );
            //Draw2D::drawLine(pi, pos[j]);
        }
        force[i].add(f);
    }
}

void add_ineraction_forces_cells( int n, Vec2f* pos, Vec2f* force, int* c2p ){
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    for(int iy=1; iy<ny_; iy++){
        for(int ix=1; ix<nx_; ix++){
            int i      = nx*iy+ix;
            int iStart = c2p[i  ];
            int iEnd   = c2p[i+1];
            if( iStart==iEnd ) continue;
            int j;
            // upper row (iy-1)

            //DEBUG_counter = 0;
            //DEBUG_counter += iEnd-iStart;

            j=i-nx-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i   -1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
                      interact_cell( iStart, iEnd,                   pos, force );
            j=i   +1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i+nx-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i+nx  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i+nx+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
        }
    }
    //exit(0);
}

void clean_array( int n, Vec2f * arr ){ for(int i=0; i<n; i++){ arr[i].set(0.0); } }

double checkDiff( int n, Vec2f * ps, Vec2f * p0s ){
    double errSum2 = 0;
    for(int i=0; i<n; i++){
        Vec2f d; d.set_sub( ps[i], p0s[i] );
        double err2 = d.norm2();
        errSum2    += err2;
        if ( err2 > 1e-8 ){
            printf( "%i is (%g,%g,) should be (%g,%g,) \n", i, ps[i].x, ps[i].y, p0s[i].x, p0s[i].y );
            exit(0);
            break;

        }
    }
    return errSum2;
}

void add_confine_force( const Vec2f& center, float strength ){
    for(int i=0; i<n; i++){ force[i].add_mul( pos[i], strength ); }
}

void add_external_force( int n, Vec2f * force, const Vec2f& center, float strength, float width ){
    float w2 = width*width;
    for(int i=0; i<n; i++){
        Vec2f dp; dp.set_sub( pos[i], center );
        //float ir2 = 1/(dp.norm2() + w2);
        //float ir  = sqrt(ir2);
        //float fr = strength*ir2*ir;
        if( dp.norm2()>w2 ) continue;
        force[i].add_mul( dp, strength );
    }
}

inline void wrap( Vec2f& p ){
    bool wrapped = false;
    if     (p.x> xspan){ p.x = -2*xspan+p.x; wrapped = true; }
    else if(p.x<-xspan){ p.x =  2*xspan+p.x; wrapped = true; }
    if     (p.y> yspan){ p.y = -2*yspan+p.y; wrapped = true; }
    else if(p.y<-yspan){ p.y =  2*yspan+p.y; wrapped = true; }
    //if(wrapped){
    //    //v.mul(0.5);
    //    double vr2 = v.norm2();
    //    if(vr2>1e+4) v.set(0.0);
    //    if(vr2>1e+2) v.mul(0.5);
    //}
    //double vr2 = v.norm2();
    //if(vr2>1) v.set(0.5);
}

float move_leap_frog( int n, Vec2f * pos, Vec2f * vel, Vec2f * force, float dt ){
    float cdamp = 1 - damp*dt;
    double f2max = 0;
    for(int i=0; i<n; i++){
        float f2 = force[i].norm2();
        if(f2>f2max) f2max = f2;
    }
    if(f2max > F2MAX ) dt *= sqrt( sqrt(F2MAX/f2max) );
    for(int i=0; i<n; i++){
        vel[i].mul(cdamp);
        vel[i].add_mul( force[i], dt );
        pos[i].add_mul( vel[i]  , dt );
        wrap( pos[i] );
    }
    return f2max;
}

#endif

