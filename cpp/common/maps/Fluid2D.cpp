
#include "Fluid2D.h"  // THE HEADER

void Fluid2D::allocate( Vec2i ns_ ){
    setN(ns_);
    //_realloc(vel  ,ntot);
    //_realloc(vel_ ,ntot);
    _realloc(vx,ntot);   _realloc(vx_,ntot);
    _realloc(vy,ntot);   _realloc(vy_,ntot);
    _realloc(dens,ntot); _realloc(dens_,ntot);
    _realloc(div ,ntot);
    _realloc(p   ,ntot);
    _realloc(source,ntot);

    setarr(ntot,vx,0);   setarr(ntot,vx_,0);
    setarr(ntot,vy,0);   setarr(ntot,vy_,0);
    setarr(ntot,dens,0); setarr(ntot,dens_,0);
    setarr(ntot,div,0);
    setarr(ntot,p  ,0);
    setarr(ntot,dens_,0);
};

 void Fluid2D::corners(double* u){
    //u[ 0 ][ 0 ] = 0.5 * (u[1][ 0 ] + u[ 0 ][1]);
    //u[ 0 ][N+1] = 0.5 * (u[1][N+1] + u[ 0 ][N]);
    //u[N+1][ 0 ] = 0.5 * (u[N][ 0 ] + u[N+1][1]);
    //u[N+1][N+1] = 0.5 * (u[N][N+1] + u[N+1][N]);

    /*
    iymax = (n.y-1)*n.x;
    u[0     ] = 0.5*( u[1] + u[n.x] );
    u[iymax ] = 0.5*( u[iymax+1] + u[iymax-n.x]);
    u[n.x-1 ] = 0.5*( u[n.x-2] + u[n.x-1] );
    u[ntot-1] = 0.5*( u[n.x-2+iymax] + u[1+iymax]);
    */
}

void Fluid2D::boundary_zero(double* u) {   // refflective boundary condition
    //int ioff = n.x*(n.y-1); for(int i=1;   i<n.x;  i++   ){ u[i]=0; u[i+ioff]=0; }
    //int ioff = n.x*(n.y-1); for(int i=n.x; i<ntot; i+=n.x){ u[i]=0; u[i+ioff]=0; }
    for(int i=1; i<n.x-1; i++ ){ u[ip2i({i,0})]=0; u[ip2i({i,n.y-1})]=0; }
    for(int i=1; i<n.y-1; i++ ){ u[ip2i({0,i})]=0; u[ip2i({n.x-1,i})]=0; }
    //corners(u);
}

// refflective boundary condition
void Fluid2D::boundary_reflect( double* u) {
    //for(int i = 1; i <= N; i++) { x[0][i] = -x[1][i]; x[N+1][i] = -x[N][i]; }
    //for(int i = 1; i <= N; i++) { x[i][0] = -x[i][1]; x[i][N+1] = -x[i][N]; }
    for(int i=1; i<n.x-1; i++ ){ u[ip2i({i,0})] = -u[ip2i({i,1})]; u[ip2i({i,n.y-1})]=-u[ip2i({i,n.y-2})]; }
    for(int i=1; i<n.y-1; i++ ){ u[ip2i({0,i})] = -u[ip2i({1,i})]; u[ip2i({n.x-1,i})]=-u[ip2i({n.y-2,i})]; }
    //corners(u);
}

// absorbing boundary condition
void Fluid2D::boundary_absorb( double* u) {
    //for(int i=1; i<=N; i++){ x[0][i]= x[1][i]; x[N+1][i] = x[N][i]; }
    //for(int i=1; i<=N; i++){ x[i][0]= x[i][1]; x[i][N+1] = x[i][N]; }
    for(int i=1; i<n.x-1; i++ ){ u[ip2i({i,0})] = u[ip2i({i,1})]; u[ip2i({i,n.y-1})]=u[ip2i({i,n.y-2})]; }
    for(int i=1; i<n.y-1; i++ ){ u[ip2i({0,i})] = u[ip2i({1,i})]; u[ip2i({n.x-1,i})]=u[ip2i({n.y-2,i})]; }
    //corners(u);
}

 //======== periodic boundary condition
void Fluid2D::boundary_periodic( double* u) {
    //for(int i = 1; i <= N; i++){  u[0][i] = u[N][i]; u[N+1][i] = u[1][i]; }
    //for(int i = 1; i <= N; i++){  x[i][0] = x[i][N]; x[i][N+1] = x[i][1]; }
    for(int i=1; i<n.x-1; i++ ){ u[ip2i({i,0})] = u[ip2i({i,n.y-2})]; u[ip2i({i,n.y-1})]=u[ip2i({i,1})]; }
    for(int i=1; i<n.y-1; i++ ){ u[ip2i({0,i})] = u[ip2i({n.x-2,i})]; u[ip2i({n.x-1,i})]=u[ip2i({1,i})]; }
    //corners(u);
}

//======== set corner points as average of line

void Fluid2D::set_bnd( int b, double* u) {
    switch(b){
        case 0: boundary_zero    (u);
        case 1: boundary_reflect (u);
        case 2: boundary_absorb  (u);
        case 3: boundary_periodic(u);
    }
    //corners(u);
}

void Fluid2D::advect(double* u, double* u0, double* vx, double* vy, double dt) {
//	int i0, j0, i1, j1;
//	double x, y, s0, t0, s1, t1;
    for (int iy=1; iy<n.y-1; iy++) {
        for (int ix=1; ix<n.x-1; ix++) {
            int i = ip2i( {ix,iy} );
            double x = ix - dt * n.x * vx[i];   // how far it flows ?
            double y = iy - dt * n.y * vy[i];
            //Draw2D::drawLine({ix*0.1,iy*0.1},{x*0.1,y*0.1});
            u[i] = interpBilinear( {x,y}, u0 );
        }
    }
}

void Fluid2D::diverg( double* div, double* p, double* vx, double* vy ){
    dx     = 1.0/(n.x-2);
    dy     = 1.0/(n.y-2);
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
        int i = ip2i({ix,iy});
        div[i] = -0.5*(                                        // compute velocity divergence
               dx*( vx[i+1  ]-vx[i-1  ] )
            +  dy*( vy[i+n.x]-vy[i-n.x] )
            );
        p[i] = 0;
    }}
    boundary_absorb(div);
    boundary_absorb(p);
}

void Fluid2D::pressureBlur( double* div, double* p ){
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
        int i = ip2i({ix,iy});
        p[i] = 0.25*(  div[i]          // evaluate pressure from velocity divergence
            + p[i-1  ] + p[i+1  ]        // average pressure, to blur
            + p[i-n.x] + p[i+n.x]
        );
    }}
    boundary_absorb(p);
}

void Fluid2D::accelerate( double* p, double* vx, double* vy, double dt ){
    double dxInv=-0.5*dt/dx;
    double dyInv=-0.5*dt/dy;
    float fsc = 1000.0;
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
        int i = ip2i({ix,iy});
        double fx = p[i+1  ]-p[i-1  ];
        double fy = p[i+n.x]-p[i-n.x];
        //if( (fx*fx+fy*fy)*fsc*fsc > 0.0001 )  Draw2D::drawLine({ix*0.1,iy*0.1},{ix*0.1+fx*fsc,iy*0.1+fy*fsc});
        vx[i] += dxInv * fx;
        vy[i] += dyInv * fy;
    }}
    boundary_reflect(vx); //  vx[:][] reflective
    boundary_reflect(vy); //  vx[:][] reflective
}

void Fluid2D::diffuse( double* u, double* u0, double diff, double dt ){    // diffuse velocity
    double a = dt * diff*ntot;
    double a_frac = 1/(1+4*a);
    //printf( " a %f a_frac %f  diff %f dt %f ntot %i \n", a, a_frac, diff, dt, ntot );
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
            int i = ip2i( {ix,iy} );
            //u[i] = ( u0[i] + ( u[i-1]+u[i+1] + u[i-n.x]+u[i+n.x] )*a )*a_frac; // we don't need index warping here
            u[i] = ( u0[i] + ( u0[i-1]+u0[i+1] + u0[i-n.x]+u0[i+n.x] )*a )*a_frac; // we don't need index warping here
    } }
}

void Fluid2D::fluidStep(double dt){

    for(int i=0;i<ndiffuse;i++){ diffuse(vx_, vx, visc, dt); SWAP(vx_,vx,double*); }
    for(int i=0;i<ndiffuse;i++){ diffuse(vy_, vy, visc, dt); SWAP(vy_,vy,double*); }
    //boundary_reflect(vx);
    //boundary_reflect(vy);
    /*
    diverg( div, p, vx, vy );
    //for(int i=0;i<npressure;i++) pressureBlur( div, p );
    for(int i=0;i<npressure;i++){ diffuse(p,div, visc, dt); SWAP(p,div,double*); }
    accelerate( p, vx, vy );
    */
    advect(vx, vx_, vx_, vy_, dt); //SWAP(vx_,vx,double*); //boundary_reflect(vx);
    advect(vy, vy_, vx_, vy_, dt); //SWAP(vy_,vy,double*); //boundary_reflect(vy);
    //acum( ntot, vx_, vx, dt );
    //acum( ntot, vy_, vy, dt );
    //SWAP(vx_,vx,double*);
    //SWAP(vy_,vy,double*);

    diverg( div, p, vx, vy ); boundary_absorb(div);
    //for(int i=0;i<npressure;i++) pressureBlur( div, p );
    for(int i=0;i<npressure;i++){ diffuse(p,div, visc, dt); SWAP(p,div,double*); }
    //SWAP(p,div,double*);
    accelerate( p, vx, vy, 1.0 );

    // ==== vel_step
    /*
    acum( ntot, vx, vx_, dt );
    acum( ntot, vy, vy_, dt );
    for(int i=0;i<ndiffuse;i++){ diffuse(vx, vx_, visc, dt); boundary_reflect(vx); }
    for(int i=0;i<ndiffuse;i++){ diffuse(vy, vy_, visc, dt); boundary_reflect(vy); }
    SWAP(vx_,vx,double*);
    SWAP(vy_,vy,double*);
    //project(vx, vy, vx0, vy0);
    diverg( div, p, vx, vy );
    for(int i=0;i<npressure;i++) pressureBlur( div, p );
    accelerate( p, vx, vy );   SWAP(vx_,vx,double*); SWAP(vy_,vy,double*);
    advect(vx, vx_, vx_, vy_, dt); boundary_reflect(vx);
    advect(vy, vy_, vx_, vy_, dt); boundary_reflect(vy);
    //project(vx, vy, vx0, vy0);             // meaning of variables:      void project(vx, vy, p, diverg
    diverg( div, p, vx, vy );
    for(int i=0;i<npressure;i++) pressureBlur( div, p );
    accelerate( p, vx, vy );
    // === dens step
    acum    (ntot,dens,dens_,dt);  SWAP(dens_,dens,double*);
    for(int i=0;i<ndiffuse;i++){ diffuse (dens ,dens_, diff, dt); SWAP(dens_, dens,double*); }
    boundary_absorb(dens);
    advect  (dens, dens_, vx, vy, dt);
    boundary_absorb(dens);
    */
}


/*
void project( Vec2d* v, double* p, double* diverg) {
    double hx     = 1.0/n.y;
    double hy     = 1.0/n.y;
    //double h_frac = 1.0/h;
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
        int i = ip2i( {ix,iy} );
        diverg[i] = -0.5*(                                        // compute velocity divergence
               hx*(v[i+1  ].x - v[i-1  ].x )
            +  hy*(v[i+n.x].y - v[i-n.x].y )
            );
        p[i] = 0;
    }}
    set_bnd(0, diverg);
    set_bnd(0, p);
    for (int k = 0; k < npressure; k++){
        for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
            int i = ip2i( {ix,iy} );
            p[i] = 0.25*(  diverg[i]          // evaluate pressure from velocity divergence
                 + p[i-1  ] + p[i+1  ]        // average pressure, to blur
                 + p[i-n.x] + p[i+n.x]
            );
        }}
        set_bnd(0, p);
    }
    acum( ntot, p, source, dt );
    for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
        int i = ip2i( {ix,iy} );
        v[i].add(
            -0.5 * n.x * (p[i+1  ] - p[i-1]  ),
            -0.5 * n.y * (p[i+n.x] - p[i-n.x])
        );
        //vx[i][j] -= 0.5 * h_frac * (p[i+1][ j ] - p[i-1][ j ]);         // accelerate mass
        //vy[i][j] -= 0.5 * h_frac * (p[ i ][j+1] - p[ i ][j-1]);
    }}
    //set_bnd(1, vx); //  vx[:][] reflective
    //set_bnd(2, vy); //  vx[:][] reflective
}
*/

//             WARP Iterpolator
/*
void interpolatedMove(double x, double y) {
x*=N;  y*=N;
int ix = (int)( x + 10000) - 10000;  // fast floor
int iy = (int)( y + 10000) - 10000;
double wx = x - ix; double mx = 1.0 - wx;
double wy = y - iy; double my = 1.0 - wy;
//if((abs(ix-(N/2))>(N/2-2))||(abs(iy-(N/2))>(N/2-2))) println(  x+" "+y+" "+ix+" "+iy);
//ix++; iy++;
dx = my* (mx*vx[ix][iy]+ wx*vx[ix+1][iy]) + wy*(mx*vx[ix][iy+1]+ wx*vx[ix+1][iy+1])  ;
dy = my* (mx*vy[ix][iy]+ wx*vy[ix+1][iy]) + wy*(mx*vy[ix][iy+1]+ wx*vy[ix+1][iy+1])  ;
}
*/
