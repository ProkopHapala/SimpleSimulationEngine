
#ifndef Grid_h
#define Grid_h

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Vec3.h"
#include "Mat3.h"
//#include <string.h>

#include <vector>

// ================= MACROS

#define fast_floor_offset  1000
#define fast_floor( x )    ( ( (int)( (x) + fast_floor_offset ) ) - fast_floor_offset )
#define i3D( ix, iy, iz )  ( (iz)*nxy + (iy)*nx + (ix)  )

// ================= CONSTANTS

// Force-Field namespace
class GridShape {
	public:
	Vec3d   pos0;
	Mat3d   cell;       // lattice vector
	Mat3d   dCell;      // basis vector of each voxel ( lattice vectors divided by number of points )
	Mat3d   diCell;     // inversion of voxel basis vector
	Vec3i   n;          // number of pixels along each basis vector

	//inline Vec3d * allocateArray_Vec3d(){ return new Vec3d[n.x*n.y*n.z); }
	inline int getNtot() const {return n.x*n.y*n.z ; }

	inline void updateCell(){
        dCell.a.set_mul( cell.a, 1.0d/n.a );
		dCell.b.set_mul( cell.b, 1.0d/n.b );
		dCell.c.set_mul( cell.c, 1.0d/n.c );
		dCell.invert_T_to( diCell );
	}

	inline void setCell( const Mat3d& cell_ ){
		//n.set( n_ );
		cell.set( cell_ );
        updateCell();
	}

	//inline void set( int * n_, double * cell_ ){ set( *(Vec3d*) n_, *(Mat3d*)cell_ ); };

	inline void grid2cartesian( const Vec3d& gpos, Vec3d& cpos ) const {
		cpos.set_mul( dCell.a, gpos.x );
		cpos.add_mul( dCell.b, gpos.y );
		cpos.add_mul( dCell.c, gpos.z );
	}

	inline void cartesian2grid( const Vec3d& cpos, Vec3d& gpos ) const {
		gpos.a = cpos.dot( diCell.a );
		gpos.b = cpos.dot( diCell.b );
		gpos.c = cpos.dot( diCell.c );
	}

	void printCell() const {
	    printf( " n      %i %i %i \n", n.x,        n.y,       n.z        );
	    printf( " a      %f %f %f \n", cell.a.x,   cell.a.y,   cell.a.z  );
	    printf( " b      %f %f %f \n", cell.b.x,   cell.b.y,   cell.b.z  );
	    printf( " c      %f %f %f \n", cell.c.x,   cell.c.y,   cell.c.z  );
	    printf( " da     %f %f %f \n", dCell.a.x,  dCell.a.y,  dCell.a.z  );
	    printf( " db     %f %f %f \n", dCell.b.x,  dCell.b.y,  dCell.b.z  );
	    printf( " dc     %f %f %f \n", dCell.c.x,  dCell.c.y,  dCell.c.z  );
	    printf( " inv_da %f %f %f \n", diCell.a.x, diCell.a.y, diCell.a.z );
	    printf( " inv_db %f %f %f \n", diCell.b.x, diCell.b.y, diCell.b.z );
	    printf( " inv_dc %f %f %f \n", diCell.c.x, diCell.c.y, diCell.c.z );
    }

    /*
    void saveXSF( char * fname, Vec3d * FF, int icomp ){
        printf( "saving %s\n", fname );
        FILE *fout;
        fout = fopen(fname,"w");
        fprintf( fout, "   ATOMS\n" );
        fprintf( fout, "    1   0.0   0.0   0.0\n" );
        fprintf( fout, "\n" );
        fprintf( fout, "BEGIN_BLOCK_DATAGRID_3D\n" );
        fprintf( fout, "   some_datagrid\n" );
        fprintf( fout, "   BEGIN_DATAGRID_3D_whatever\n" );
        fprintf( fout, "%i %i %i\n", n.x, n.y, n.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", pos0.x,   pos0.x,   pos0.x   );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        for ( int ic=0; ic<nz; ic++ ){
            for ( int ib=0; ib<ny; ib++ ){
                for ( int ia=0; ia<nx; ia++ ){
                   int i = i3D( ia, ib, ic );
                   fprintf( fout, "%6.5e\n", ((double*)(FF+i))[icomp] );
                }
            }
        }
        fprintf( fout, "   END_DATAGRID_3D\n" );
        fprintf( fout, "END_BLOCK_DATAGRID_3D\n" );
        fclose(fout);
    }
    */

    void toXSF( FILE* fout, Vec3d * FF, int icomp ) const {
        fprintf( fout, "BEGIN_BLOCK_DATAGRID_3D\n" );
        fprintf( fout, "   some_datagrid\n" );
        fprintf( fout, "   BEGIN_DATAGRID_3D_whatever\n" );
        fprintf( fout, "%i %i %i\n", n.x, n.y, n.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", pos0.x,   pos0.x,   pos0.x   );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
        fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
        int nx  = n.x; 	int ny  = n.y; 	int nz  = n.z; int nxy = ny * nx;
        for ( int ic=0; ic<nz; ic++ ){
            for ( int ib=0; ib<ny; ib++ ){
                for ( int ia=0; ia<nx; ia++ ){
                   int i = i3D( ia, ib, ic );
                   fprintf( fout, "%6.5e\n", ((double*)(FF+i))[icomp] );
                }
            }
        }
        fprintf( fout, "   END_DATAGRID_3D\n" );
        fprintf( fout, "END_BLOCK_DATAGRID_3D\n" );
    }

    void saveXSF( char * fname, Vec3d * FF, int icomp )const {
        printf( "saving %s\n", fname );
        FILE *fout;
        fout = fopen(fname,"w");
        fprintf( fout, "   ATOMS\n" );
        fprintf( fout, "    1   0.0   0.0   0.0\n" );
        fprintf( fout, "\n" );
        toXSF( fout, FF, icomp );
        fclose(fout);
    }

};

// interpolation of vector force-field Vec3d[ix,iy,iz] in periodic boundary condition
inline double interpolate3DWrap( double * grid, const Vec3i& n, const Vec3d& r ){
	int xoff = n.x<<3; int imx = r.x +xoff;	double tx = r.x - imx +xoff;	double mx = 1 - tx;		int itx = (imx+1)%n.x;  imx=imx%n.x;
	int yoff = n.y<<3; int imy = r.y +yoff;	double ty = r.y - imy +yoff;	double my = 1 - ty;		int ity = (imy+1)%n.y;  imy=imy%n.y;
	int zoff = n.z<<3; int imz = r.z +zoff;	double tz = r.z - imz +zoff;	double mz = 1 - tz;		int itz = (imz+1)%n.z;  imz=imz%n.z;
	int nxy = n.x * n.y; int nx = n.x;
	//double out = grid[ i3D( imx, imy, imz ) ];

	//double out = mz * my * (  ( mx * grid[ i3D( imx, imy, imz ) ] ) +  ( tx * grid[ i3D( itx, imy, imz ) ] );
	//ty * ( mx * grid[ i3D( imx, ity, imz ) ] ) +  ( tx * grid[ i3D( itx, ity, imz ) ] ) );

	double out = mz * (
	my * ( ( mx * grid[ i3D( imx, imy, imz ) ] ) +  ( tx * grid[ i3D( itx, imy, imz ) ] ) ) +
	ty * ( ( mx * grid[ i3D( imx, ity, imz ) ] ) +  ( tx * grid[ i3D( itx, ity, imz ) ] ) ) )
               + tz * (
	my * ( ( mx * grid[ i3D( imx, imy, itz ) ] ) +  ( tx * grid[ i3D( itx, imy, itz ) ] ) ) +
	ty * ( ( mx * grid[ i3D( imx, ity, itz ) ] ) +  ( tx * grid[ i3D( itx, ity, itz ) ] ) ) );
	return out;
}

// interpolation of vector force-field Vec3d[ix,iy,iz] in periodic boundary condition
inline Vec3d interpolate3DvecWrap( Vec3d * grid, const Vec3i& n, const Vec3d& r ){
	int xoff = n.x<<3; int imx = r.x +xoff;	double tx = r.x - imx +xoff;	double mx = 1 - tx;		int itx = (imx+1)%n.x;  imx=imx%n.x;
	int yoff = n.y<<3; int imy = r.y +yoff;	double ty = r.y - imy +yoff;	double my = 1 - ty;		int ity = (imy+1)%n.y;  imy=imy%n.y;
	int zoff = n.z<<3; int imz = r.z +zoff;	double tz = r.z - imz +zoff;	double mz = 1 - tz;		int itz = (imz+1)%n.z;  imz=imz%n.z;
	int nxy = n.x * n.y; int nx = n.x;
	//printf( " %f %f %f   %i %i %i \n", r.x, r.y, r.z, imx, imy, imz );
	double mymx = my*mx; double mytx = my*tx; double tymx = ty*mx; double tytx = ty*tx;
	Vec3d out;
	out.set_mul( grid[ i3D( imx, imy, imz ) ], mz*mymx );   out.add_mul( grid[ i3D( itx, imy, imz ) ], mz*mytx );
	out.add_mul( grid[ i3D( imx, ity, imz ) ], mz*tymx );   out.add_mul( grid[ i3D( itx, ity, imz ) ], mz*tytx );
	out.add_mul( grid[ i3D( imx, ity, itz ) ], tz*tymx );   out.add_mul( grid[ i3D( itx, ity, itz ) ], tz*tytx );
	out.add_mul( grid[ i3D( imx, imy, itz ) ], tz*mymx );   out.add_mul( grid[ i3D( itx, imy, itz ) ], tz*mytx );
	return out;
}

// iterate over field
//template< void FUNC( int ibuff, const Vec3d& pos_, void * args ) >
template<typename FUNC>
void interateGrid3D( const GridShape& grid, FUNC func ){
	int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z;
	//int nx  = n.z; 	int ny  = n.y; 	int nz  = n.x;
	int nxy = ny * nx;
	printf( "interateGrid3D nx,y,z (%i,%i,%i) nxy %i\n", nx,ny,nz, nxy );
	Vec3d pos;  pos.set( grid.pos0 );
	//printf(" interateGrid3D : args %i \n", args );
	for ( int ic=0; ic<nz; ic++ ){
        std::cout << "ic " << ic;
        std::cout.flush();
        std::cout << '\r';
		for ( int ib=0; ib<ny; ib++ ){
	        for ( int ia=0; ia<nx; ia++ ){
			    int ibuff = i3D( ia, ib, ic );
                //FUNC( ibuff, {ia,ib,ic}, pos );
                //pos = pos0 + dCell.c*ic + dCell.b*ib + dCell.a*ia;
                func( ibuff, pos );
                //printf("(%i,%i,%i)(%3.3f,%3.3f,%3.3f)\n",ia,ib,ic,pos.x,pos.y,pos.z);
				pos.add( grid.dCell.a );
			}
			pos.add_mul( grid.dCell.a, -nx );
			pos.add( grid.dCell.b );
		}
		//exit(0);
		pos.add_mul( grid.dCell.b, -ny );
		pos.add( grid.dCell.c );
	}
    printf ("\n");
}

void getIsovalPoints_a( const GridShape& grid, double isoval, Vec3d  *FF, std::vector<Vec3d>& points ){
    int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z; int nxy = ny * nx;
    int ii = 0;
    for ( int ic=0; ic<(nz-1); ic++ ){
        for ( int ib=0; ib<ny; ib++ ){
            for ( int ia=0; ia<nx; ia++ ){
                int     i1 = i3D( ia, ib,  ic    );
                int     i2 = i3D( ia, ib, (ic+1) );
                double df1 = FF[i1].z-isoval;
                double df2 = FF[i2].z-isoval;
                //if(ii<1000)printf( " %i (%i,%i,%i) (%g,%g)\n", ii, ia,ib,ic, df1, df2 );
                if( (df1*df2)<0 ){
                    double fc = df1/(df1-df2);
                    points.push_back( grid.dCell.a*ia + grid.dCell.b*ib + grid.dCell.c*(ic+fc) );

                    int ip = points.size()-1;
                    if( ip < 1000 ) printf( " %i (%i,%i,%i) (%g,%g,%g)\n", ip, ia,ib,ic, points[ip].x, points[ip].y, points[ip].z );
                }
                ii++;
            }
        }
    }
}

void getIsoSurfZ( const GridShape& grid, double isoval, bool sign, Vec3d  *FF, Vec3d *pos, Vec3d * normal ){
    int nx  = grid.n.x; 	int ny  = grid.n.y; 	int nz  = grid.n.z; int nxy = ny * nx;
    int ii = 0;
    //printf("%i %i %i \n", nx,ny,nxy );
    for ( int ib=0; ib<ny; ib++ ){
        for ( int ia=0; ia<nx; ia++ ){
            int ibuff = i3D( ia, ib,  0 );
            double ofz = FF[ibuff].x;
            double fz;
            int ic;
            //printf( "iba (%i,%i)\n", ib,ia );
            for ( ic=nz-1; ic>1; ic-- ){
                int ibuff_ = ibuff + nxy*ic;
                fz = FF[ibuff_].z;
                if( (fz>isoval)==sign ){
                    ibuff = ibuff_;
                    break;
                }
                ofz = fz;
            }
            //double fc = (ofz-isoval)/(ofz-fz);
            double fc = 1-((ofz-isoval)/(ofz-fz));
            //double fc = 0;
            int ibxy  = ib*nx + ia;
            //printf( "ibxy %i %i \n", ibxy, ibuff );
            pos   [ibxy] = grid.dCell.a*ia + grid.dCell.b*ib + grid.dCell.c*(ic+fc);
            //DEBUG
            //normal[ibxy] = FF[ibuff-1];
            normal[ibxy] = FF[ibuff];
            //DEBUG
            normal[ibxy].normalize();
            //normal[ibxy] = interpolate3DvecWrap( FF, grid.n, {ia,ib, ic+fc } );
            //printf("done \n");
        }
    }
}

void writePrimCoord( FILE* fout, Mat3d& cell, int natoms, Vec3d * apos, int * iZs ){
    fprintf( fout, "CRYSTAL\n" );
    fprintf( fout, "PRIMVEC\n" );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
    fprintf( fout, "CONVVEC\n" );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.a.x, cell.a.y, cell.a.z );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.b.x, cell.b.y, cell.b.z );
    fprintf( fout, "%5.10f %5.10f %5.10f \n", cell.c.x, cell.c.y, cell.c.z );
    fprintf( fout, "PRIMCOORD\n" );
    fprintf( fout, "%i %i\n", natoms, 1 );
    for(int i=0; i<natoms; i++){
        fprintf( fout, "%i %3.8f %3.8f %3.8f\n", iZs[i], apos[i].x, apos[i].y, apos[i].z );
    }
}

void saveXSF( char * fname, GridShape& grid, Vec3d * FF, int icomp, int natoms, Vec3d * apos, int * iZs ){
    printf( "saving %s\n", fname );
    FILE *fout;
    fout = fopen(fname,"w");
    writePrimCoord( fout, grid.cell, natoms, apos, iZs );
    grid.toXSF    ( fout, FF, icomp );
    fclose(fout);
}


#endif









