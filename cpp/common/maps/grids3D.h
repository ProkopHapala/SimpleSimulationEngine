#ifndef  grids3D_h
#define  grids3D_h

/*

This is header-only module with mostly function templates for common operations on regular grids in 3D

we try to avoid making complex classes in order to keep things simple, modular and independent ( with classes we usually and up in dependency-hell and too complicated templates )  

operations such as:
  - support rutines for particle in cell (PIC) sumulation 
  - nearest neighbor interaction in 3D accelerated by grid
  - boundary condition resolution for 3D grid

*/

#include "integerOps.h"
#include "fastmath.h"
#include "Vec3.h"

class GridRulerInterface{
    public:
};

class CubeGridRuler : public GridRulerInterface {
    public:
    double step;
    double invStep;
    Vec3d  pos0;
    
    inline setStep( double step_ ){ step=step_; invStep=1/step; };
    
    inline void pos2box( const Vec3d& pos, Vec3i& ipos, Vec3d& dpos ) const {
        dpos.x = x2grid( pos.x-pos0.x, step, invStep, ipos.x );
        dpos.y = x2grid( pos.y-pos0.y, step, invStep, ipos.y );
        dpos.z = x2grid( pos.z-pos0.z, step, invStep, ipos.z );
        //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", pos.x, pos.y, pos.z, ipos.x,ipos.y,ipos.z);
    }
};


namespace Grids3D {

    // should we move this to Vec3.h ?
    inline wrapIndex( Vec3i& ipos, const Vec3i& n ){
        if(ipos>0){}else{};
    }


    //typedef KeyType  unit32_t;
    //typedef ValType  unit32_t;
    
    //   search 2x2x2 neighborhood of 
    template< OBJECT o, void insert_func( int ix, int iy, int iz, int val) >
    void insert_SphereOfInfluence( OBJECT o, const Vec3i& ipos, const Vec3d& dpos, double r ){
	    insert( o, ipos.x, ipos.y, ipos.z );
	    int    dix =0,diy =0,diz =0;
	    double dr2x=0,dr2y=0,dr2z=0;
	    double mr = 1-r;
	    if     (  dpos.x < r  ){ insert_func( ipos.x-1, ipos.y  , ipos.z, o );    dix=-1; dr2x = sq(  dpos.x); }
	    else if(  dpos.x > mr ){ insert_func( ipos.x+1, ipos.y  , ipos.z, o );    dix=+1; dr2x = sq(1-dpos.x); }
	    if     (  dpos.y < r  ){ insert_func( ipos.x  , ipos.y-1, ipos.z, o );    diy=-1; dr2y = sq(  dpos.y); }
	    else if(  dpos.y > mr ){ insert_func( ipos.x  , ipos.y+1, ipos.z, o );    diy=+1; dr2y = sq(1-dpos.y); }
	    if     (  dpos.z < r  ){ insert_func( ipos.x  , ipos.y  , ipos.z-1, o );  diz=-1; dr2z = sq(  dpos.z); }
	    else if(  dpos.z > mr ){ insert_func( ipos.x  , ipos.y  , ipos.z+1, o );  diz=+1; dr2z = sq(1-dpos.z); }

	    double r2 = r*r;
	    if ( dr2x+dr2y      < r2 ){ insert_func( ipos.x+dix, ipos.y+diy, ipos.z o      ); }
	    if ( dr2x+dr2z      < r2 ){ insert_func( ipos.x+dix, ipos.y    , ipos.z+diz, o ); }
	    if ( dr2y+dr2z      < r2 ){ insert_func( ipos.x    , ipos.y+diy, ipos.z+diz, o ); }
	    if ( dr2x+dr2y+dr2z < r2 ){ insert_func( ipos.x+dix, ipos.y+diy, ipos.z+diz, o ); }

	    //if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
	    //printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
    }
    
    template< OBJECT o, void insert_func( int ix, int iy, int iz, int val) >
    inline void insert( OBJECT o, const CubeGridRuler& ruler, const Vec3d& pos, double r ){
        Vec3d dpos; Vec3i ipos;
        pos2box( pos, ipos,dpos );
        insert<insert_func>( val, ipos, dpos, r );
    }

}

#endif
