#ifndef  NBody_h
#define  NBody_h

#include "OCL.h"

#define R2SAFE  1.0e-2
#define F2MAX   10.0

constexpr float R_MAX  = 1.8;
constexpr float R2MAX = R_MAX*R_MAX;

//int n = 1024;
//constexpr int n = 32;
constexpr int n = 128;
float time_step = 0.05;
//constexpr int n = 1024*2;
float damp      = 0.5;
Vec2f pos   [n];
Vec2f pos_  [n];
Vec2f vel   [n];
Vec2f force [n];
Vec2f force_[n];


// cell acceleration buffer
int constexpr nx    = 16;
int constexpr ny    = 16;
int constexpr ncell = nx*ny;
int           cell2pos  [ncell+1];
int           cellNs    [ncell];
//int           cellNcount[ncell];
int           atom2cell[n];
float cell_size     = 2.0;
float grid_orig [2] = {-16.0,-16.0};



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
    float ir2 = 1/( dp.norm2() + R2SAFE );
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
                pos[i].set( ix*step -0.5*nx + randf(-drnd,drnd), iy*step -0.5*ny +randf(-drnd,drnd) );
                vel  [i].set(0.0);
                force[i].set(0.0);
            }
            i++;
        }
    }
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
    /*
    for(int i=0; i<n; i++){
        Vec2f& p  = pos_[i];
        int ix    = (p.x - grid_orig[0])*inv_size;
        int iy    = (p.y - grid_orig[1])*inv_size;
        int icell = nx*iy + ix;
        printf( "%i (%3.3f,%3.3f) (%i,%i) %i \n", i, p.x,p.y, ix, iy, icell );
    }
    */
    //exit(0);
    for(int i=0; i<n; i++){ pos[i]=pos_[i]; }
}


void add_ineraction_forces(){
    for(int i=0; i<n; i++){
        Vec2f f; f.set(0.0);
        Vec2f& pi = pos[i];
        for(int j=0; j<n; j++){ acum_force( pi, pos[j], f ); }
        force[i] = f;
    }
}

void interact_cell( int iStart, int iEnd, int jStart, int jEnd, Vec2f* pos, Vec2f* force ){
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

void add_ineraction_forces_cells( Vec2f* pos, int* c2p ){
    const int nx_ = nx-1;
    const int ny_ = ny-1;
    for(int iy=1; iy<ny_; iy++){
        for(int ix=1; ix<nx_; ix++){
            int i      = nx*iy+ix;
            int iStart = c2p[i  ];
            int iEnd   = c2p[i+1];
            int j;
            // upper row (iy-1)

            j=i-nx-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i   -1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
                      interact_cell( iStart, iEnd,                   pos, force );
            j=i   +1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );
            j=i-nx+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );

            //printf( " (%i,%i) %i (%i,%i) (%i,%i) \n", ix, iy, i, iStart, iEnd,  c2p[j], c2p[j+1] );
/*
            if(iy>1  ){
                int i_ = i - nx;
                if(ix>1  ){ j=i_-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force ); }  //   [ iy-1 , ix-1 ]
                            j=i_  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );    //   [ iy-1 , ix   ]
                if(ix<nx_){ j=i_+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force ); }  //   [ iy-1 , ix+1 ]
            }

            // middle row (iy)
                if(ix>1  ){ j=i-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );  }  //  [ iy   , ix-1 ]
                                   interact_cell( iStart, iEnd, iStart, iEnd, pos, force );         //  [ iy   , ix   ]
                if(ix<nx_){ j=i+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );  }  //  [ iy   , ix+1 ]
            // upper row (iy+1)
            if(iy<ny_){
                int i_ = i + nx;
                if(ix>1  ){ j=i_+1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force ); }  //  [iy+1  , ix-1 ]
                            j=i_  ; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force );    //  [iy+1  , ix   ]
                if(ix<nx_){ j=i_-1; interact_cell( iStart, iEnd, c2p[j], c2p[j+1], pos, force ); }  //  [iy+1  , ix+1 ]
            }
*/
        }
    }
    //exit(0);
}

void clean_force( ){ for(int i=0; i<n; i++){ force[i].set(0.0); } }

void add_confine_force( const Vec2f& center, float strength ){
    for(int i=0; i<n; i++){ force[i].add_mul( pos[i], strength ); }
}

void add_external_force( const Vec2f& center, float strength, float width ){
    float w2 = width*width;
    for(int i=0; i<n; i++){
        Vec2f dp; dp.set_sub( pos[i], center );
        float ir2 = 1/(dp.norm2() + w2);
        float ir  = sqrt(ir2);
        force[i].add_mul( dp, strength*ir2*ir );
    }
}

#define xspan 12.0
#define yspan 12.0

inline void wrap( Vec2f& p, Vec2f& v ){
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
    double vr2 = v.norm2();
    if(vr2>1) v.set(0.5);

}

float move_leap_frog( float dt ){
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
        wrap( pos[i], vel[i] );
    }
    return f2max;
}


#endif

