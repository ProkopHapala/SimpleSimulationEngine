#ifndef TerrainHydraulics_h
#define TerrainHydraulics_h

#include "Noise.h"

class TerrainHydraulics{
public:

	int    nx, ny, ntot;
	double * ground;
	double * water;

	// countor step algorithm
	//int    nContourMax;
	int    nContour=0,nContour_old=0;
	bool   * known;
	int    * contour1;
	int    * contour2;

	// ==== inline functions

	inline int xy2i( int ix, int iy  ){ return iy<<nx + ix; };
	//inline i2xy( int i, int ix, int iy ){     };

	void allocate( int nx_, int ny_ ){
		nx = nx_; ny = ny_; ntot = nx * ny;
		ground   = new double[ntot];
		water    = new double[ntot];
        known    = new bool  [ntot]; // FIXME: this is just temporary array
	    contour1 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot
	    contour2 = new int   [ntot]; // FIXME: just temporary, 1D => size < ntot
	}

    void genTerrainNoise( int n, double scale,  double hscale,  double fdown, double strength, int seed, const Vec2d& pos0 ){
        Vec2d pos,dpos,rot,a,b;
        rot.fromAngle( seed*0.1 );
        a.set( 1.0d, 0.0d           ); a.mul(scale);
        b.set( 0.5d, 0.86602540378d ); b.mul(scale);
        double vmin=+1e+300,vmax=-1e+300;
        int ii = 0;
        for (int iy=0; iy<ny; iy++){
            for (int ix=0; ix<nx; ix++){
                pos.set( ix*a.x+iy*b.x +pos0.x, ix*a.y+iy*b.y+pos0.y );
                Noise::warpNoise3R( pos, rot, fdown, strength, n, dpos );
                //ground[ii] = 0.5*sin( 1000*dpos.x * dpos.y )+0.5;
                double val = fabs(dpos.x * dpos.y);
                ground[ii] = val;
                if(val<vmin){vmin=val;}
                if(val>vmax){vmax=val;}
                ii++;
            }
        }
        printf( " vmin %e vmax %e \n", vmin, vmax );
        double renorm = hscale/(vmax-vmin);
        for(int i=0; i<ntot; i++){ ground[i] = renorm*(ground[i]-vmin); }
    }

void init_outflow(){
    for (int i=0; i<ntot; i++){
        water[i] = 1e+300;
        known[i] = false;
    }
}

void outflow_step(){
    // clean known
    // swap old and new endpoints
    nContour_old = nContour;  nContour = 0;
    int * tmp = contour1; contour1 = contour2; contour2 = tmp;
    for ( int ii=0; ii<nContour_old; ii++ ){ known[contour1[ii]] = false; }
    // check all countor point neighbors for possible path extension
    printf( "nContour_old %i nContour %i \n", nContour_old, nContour );
    for ( int ii=0; ii<nContour_old; ii++ ){
        int    i   = contour1[ii];
        double val = water[i];
        //val +=  dval/( val - ground[i] + 0.1 ); // this is just to simulate non-zero viscosity of watter
        int iy = i/nx;
        int ix = i-(iy*nx);
        //printf( " %i %i (%i,%i)\n", ii, i, ix, iy );
        // extend in four directions, check boundary overflow
        if( ix>0      ){  extend_path( val, i, i - 1  ); }
        if( ix<(nx-1) ){  extend_path( val, i, i + 1  ); }
        if( iy>0      ){  extend_path( val, i, i - nx ); }
        if( iy<(ny-1) ){  extend_path( val, i, i + nx ); }
    }
}

void extend_path( float val, int oi, int i ){
    // evaluate objective function for proposed path
    //val += dval*ground[i];
    //printf( " %i val %f g %f w %f k %i \n", i, val, ground[i], water[i], known[i] );
    val = fmax( val, ground[i] );
    // check if new path is advantegeneous
    if( val < water[i] ){
        water[i] = val;       //  extend path
        //printf( "mod %i  \n", i );
        //gofrom[i] = oi;     //  ... optionaly we may store oi to reconstruct path backwards later
        if( !known[i] ){
                //printf( "added %i to contour \n", i );
                known[i]=true; contour2[nContour]=i; nContour++;
        } // add to list only if not already known
    }
}





};

#endif
