#ifndef Fluid2D_h
#define Fluid2D_h

#include <vector>
#include "Noise.h"
#include "arrayAlgs.h"
#include "GridIndex2D.h"
//#include "Grid2D.h"
#include "Grid2DAlgs.h"

//#include "Fluid2D.cpp"

/*
References:
http://www.dgp.toronto.edu/people/stam/reality/index.html
http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
/home/prokop/Dropbox/MyDevSW/Processing/AlgorithmicArtGeneration/PlasticOrogennyParticles
*/

void acum( int n, double* a, double* b, double dt){ for(int i=0;i<n;i++){  a[i] += b[i] * dt; } };
void swap_ptr(void*& a, void*& b){ void*& tmp=a; a=b; b=tmp; };

class Fluid2D : public Grid2DAlg { public:

    //int ndiffuse  = 20;
    //int npressure = 20;
    int ndiffuse  = 5;
    int npressure = 5;
	double N_INVERSE = (double)1.0/(double)ntot;
    double dx,dy;
    double dt;
	//double * vx     = new double[N+2][N+2];
	//double * vx_    = new double[N+2][N+2];
	//double * vy     = new double[N+2][N+2];
	//double * vy_    = new double[N+2][N+2];
    //double * dens   = new double[N+2][N+2];
    //double * dens_  = new double[N+2][N+2];

	//Vec2i ns;  // GridIndex2D
	//int ntot;  // GridIndex2D


	// QUESTION :   would it be easier to use double * vx,vy instead of Vec2d* v ????
	Vec2d  * vel    = 0;
	Vec2d  * vel_   = 0;
	double * dens   = 0;
	double * dens_  = 0;
    double * source = 0;

    // =================== Functions

    void allocate( Vec2i ns_ ){
        setN(ns_);
        _realloc(vel  ,ntot);
        _realloc(vel_ ,ntot);
        _realloc(dens ,ntot);
        _realloc(dens_,ntot);
        _realloc(source,ntot);
    };

 //---------------------------------------------------------------
 // ---------------------- Interface  UTILS-----------------------
 //---------------------------------------------------------------

 /*
	public double getDx(int x, int y) {  return vx[x+1][y+1];   }
	public double getDy(int x, int y) {  return vy[x+1][y+1]; 	}
	public void applyForce(int cellX, int cellY, double fx, double fy) {
	  vx[cellX][cellY] += fx * dt;
	  vy[cellX][cellY] += fy * dt;
	}
*/
 //---------------------------------------------------------------
 // ---------------------- MAIN ALGORITHM ------------------------
  //---------------------------------------------------------------

  void boundary_zero( double* u) {   // refflective boundary condition
  //int ioff = n.x*(n.y-1); for(int i=1;   i<n.x;  i++   ){ u[i]=0; u[i+ioff]=0; }
  //int ioff = n.x*(n.y-1); for(int i=n.x; i<ntot; i+=n.x){ u[i]=0; u[i+ioff]=0; }
  for(int i=1; i<n.x-1; i++ ){ u[ip2i({i,0})]=0; u[ip2i({i,n.y-1})]=0; }
  for(int i=1; i<n.y-1; i++ ){ u[ip2i({0,i})]=0; u[ip2i({n.x-1,i})]=0; }
}

// refflective boundary condition
void boundary_reflect( double* u) {
  //for(int i = 1; i <= N; i++) { x[0][i] = -x[1][i]; x[N+1][i] = -x[N][i]; }
  //for(int i = 1; i <= N; i++) { x[i][0] = -x[i][1]; x[i][N+1] = -x[i][N]; }
}

// absorbing boundary condition
void boundary_absorb( double* u) {
  //for(int i=1; i<=N; i++){ x[0][i]= x[1][i]; x[N+1][i] = x[N][i]; }
  //for(int i=1; i<=N; i++){ x[i][0]= x[i][1]; x[i][N+1] = x[i][N]; }
}

 //======== periodic boundary condition
void boundary_periodic( double* u) {
  //for(int i = 1; i <= N; i++){  u[0][i] = u[N][i]; u[N+1][i] = u[1][i]; }
  //for(int i = 1; i <= N; i++){  x[i][0] = x[i][N]; x[i][N+1] = x[i][1]; }
}

//======== set corner points as average of line

void corners(double* u){
/*
  u[ 0 ][ 0 ] = 0.5 * (u[1][ 0 ] + u[ 0 ][1]);
  u[ 0 ][N+1] = 0.5 * (u[1][N+1] + u[ 0 ][N]);
  u[N+1][ 0 ] = 0.5 * (u[N][ 0 ] + u[N+1][1]);
  u[N+1][N+1] = 0.5 * (u[N][N+1] + u[N+1][N]);
  */
}

void set_bnd( int b, double* u) {
 if (b==1) { boundary_reflect( u ); }
 else      { boundary_absorb(  u ); }
 //boundary_periodic( x  );
 corners(u);
}

  void diffuseVec2( Vec2d* v, Vec2d* v0, double diff) {    // diffuse velocity
		double a = dt * diff * ntot;
        double a_frac = 1/(1+4*a);
		for (int k = 0; k<ndiffuse; k++){
			for (int iy = 1; iy<n.y-1; iy++){
				for (int ix = 1; ix<n.x-1; ix++){
                    int i = ip2i( {ix,iy} );
					v[i] = ( v0[i] + ( v[i-1] + v[i+1]+ v[i-n.x]+ v[i+n.x] )*a )*a_frac; // we don't need index warping here
				}
			}
			// boundary copy ?
		}
	}

	void advect(double* d, double* d0, Vec2d* v) {
	//	int i0, j0, i1, j1;
	//	double x, y, s0, t0, s1, t1;
		for (int iy=1; iy<n.y-1; iy++) {
			for (int ix=1; ix<n.x-1; ix++) {
                int i = ip2i( {ix,iy} );
				double x = ix - dt * n.x * v[i].x;   // how far it flows ?
				double y = iy - dt * n.y * v[i].y;
                //if (x <     0.5) x =     0.5; // boundary
                //if (x > n.x + 0.5) x = n.y + 0.5;
                //if (y <     0.5) y =     0.5;
                //if (y > n.y + 0.5) y = n.y + 0.5;
                x=_clamp(x, 0.5, n.x-1.5);  // boundary condition
                y=_clamp(y, 0.5, n.y-1.5);
				int ix0 = (int)x;     //int ix1 = ix+1;
				int iy0 = (int)y;     //int iy1 = iy+1;
				int i0  = ip2i( {ix0,iy0} );
				int i1  = i0 + n.x;
				// bilinear interpolation of density
				double fx = x-ix0;  double mx = 1-fx;
				double fy = y-iy0;  double my = 1-fy;
				d[i] = my * (mx * d0[i0] + fx * d0[i0+1])
                     + fy * (mx * d0[i1] + fx * d0[i1+1]);
			}
		}
	}

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
		for (int k = 0; k < npressure; k++) {
            for (int iy=1; iy<n.y-1; iy++){ for (int ix=1; ix<n.x-1; ix++){
                int i = ip2i( {ix,iy} );
                p[i] = 0.25*(  diverg[i]                             // evaluate pressure from velocity divergence
                     + p[i-1  ] + p[i+1  ]                     // average pressure, to blur
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
// ========== Boundary conditions ===============



	void vel_step( Vec2d* v, Vec2d* v0, double visc) {
        /*
		add(vx, vx0, dt);                       add(vy, vy0, dt);
		acum( );
		SWAP(vx0, vx);                          SWAP(vy0, vy);
		diffuse(vx, vx0, visc, dt);  set_bnd(1, vx);
		diffuse(vy, vy0, visc, dt);	 set_bnd(2, vy);
		project(vx, vy, vx0, vy0);
		SWAP(vx0, vx); SWAP(vy0, vy);
		advect(vx, vx0,     vx0, vy0, dt);  set_bnd(1, vx);  // unaseni rychlosti? coze?
		advect(vy, vy0,     vx0, vy0, dt);  set_bnd(2, vy);
		project(vx, vy, vx0, vy0);             // meaning of variables:      void project(vx, vy, p, diverg)
		*/
	}

  	void dens_step(double* dens, double* dens0, double* v, double diff) {
		/*
		//add (dens , dens0 , dt);
		//SWAP       (dens0, dens  );
		//diffuse    (dens  , dens0, diff, dt);
		set_bnd(0, dens);
		SWAP       (dens0, dens  );
		advect     ( dens, dens0, vx, vy, dt);
        set_bnd(0, dens);
        */
        //acum( dens, dens0, dt );
       // advect     ( dens, dens0, vel, dt);
       // swap_ptr();
	}

	void tick(double dt, double visc, double diff) {
                //add(dens, source, dt);
//		vel_step(vel, vel_, visc, dt);              // all dynamics is solved here
//		dens_step(dens, dens_, vel, diff, dt);        // this just move mass acording to velocity
	}

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

};

#endif
