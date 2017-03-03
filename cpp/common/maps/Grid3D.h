
#ifndef  Grid3D_h
#define  Grid3D_h

//template <,N>
//class BufferTile
#include <vector>

#include "fastmath.h"
#include "Vec3.h"

template<class OBJECT>
class TILE_vector{
    std::vector<OBJECT> buff;
    inline void insert(OBJECT o){ buff.push_back(o); };
    inline void erase (        ){ buff.erase(); }
};

template<class TILE, class OBJECT, int NX_,int NY_,int NZ_>
class Grid3D{
    public:
    static constexpr int NX   = NX_;
	static constexpr int NY   = NY_;
	static constexpr int NZ   = NZ_;
	static constexpr int NXY  = NX*NY;
	static constexpr int NXYZ = NX*NY*NZ;

	Vec3d  pos0;
	double step,invStep;
	//int    count;
	//int    counts [NXYZ];
	TILE   buff   [NXYZ];

	inline int  xyz2i(        int  ix, int  iy, int  iz ){ return ix + NX*iy + NXY*iz; }
	inline void i2xyz( int i, int& ix, int& iy, int& iz ){ iz=i/NXY; i=i%NXY; iy=i/NX; ix=i%NX; } // see performance http://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder

	TILE* address( int  ix, int  iy, int  iz ){ return &buff[xyz2i(ix,iy,iz)]; };
	TILE* address( int  ixyz                 ){ return &buff[ixyz]; };

    inline void pos2box( const Vec3d& pos, Vec3i& ipos, Vec3d& dpos ){
        dpos.x = x2grid( pos.x-pos0.x, step, invStep, ipos.x );
        dpos.y = x2grid( pos.y-pos0.y, step, invStep, ipos.y );
        dpos.z = x2grid( pos.z-pos0.z, step, invStep, ipos.z );
        //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", pos.x, pos.y, pos.z, ipos.x,ipos.y,ipos.z);
    }

    inline bool validIndex( const Vec3i& ipos ){
        return (ipos.x>=0)&&(ipos.y>=0)&&(ipos.z>=0)&&(ipos.x<NX)&&(ipos.y<NY)&&(ipos.z<NZ);
    }

    // ================= Neads OBJECT

    void insert( OBJECT o, const Vec3i& ipos, const Vec3d& dpos, double r ){
		insert( o, ipos.x, ipos.y, ipos.z );
		int dix=0,diy=0,diz=0;
		double dr2x=0,dr2y=0,dr2z=0;
		double mr = 1-r;
		if     (  dpos.x < r  ){ buff[xyz2i(ipos.x-1, ipos.y  , ipos.z  )].insert( o );  dix=-1; dr2x = sq(  dpos.x); }
		else if(  dpos.x > mr ){ buff[xyz2i(ipos.x+1, ipos.y  , ipos.z  )].insert( o );  dix=+1; dr2x = sq(1-dpos.x); }
		if     (  dpos.y < r  ){ buff[xyz2i(ipos.x  , ipos.y-1, ipos.z  )].insert( o );  diy=-1; dr2y = sq(  dpos.y); }
		else if(  dpos.y > mr ){ buff[xyz2i(ipos.x  , ipos.y+1, ipos.z  )].insert( o );  diy=+1; dr2y = sq(1-dpos.y); }
		if     (  dpos.z < r  ){ buff[xyz2i(ipos.x  , ipos.y  , ipos.z-1)].insert( o );  diz=-1; dr2z = sq(  dpos.z); }
		else if(  dpos.z > mr ){ buff[xyz2i(ipos.x  , ipos.y  , ipos.z+1)].insert( o );  diz=+1; dr2z = sq(1-dpos.z); }

		double r2 = r*r;
		if ( dr2x+dr2y      < r2 ){ buff[xyz2i(ipos.x+dix, ipos.y+diy, ipos.z     )].insert( o ); }
		if ( dr2x+dr2z      < r2 ){ buff[xyz2i(ipos.x+dix, ipos.y    , ipos.z+diz )].insert( o ); }
		if ( dr2y+dr2z      < r2 ){ buff[xyz2i(ipos.x    , ipos.y+diy, ipos.z+diz )].insert( o ); }
		if ( dr2x+dr2y+dr2z < r2 ){ buff[xyz2i(ipos.x+dix, ipos.y+diy, ipos.z+diz )].insert( o ); }

		//if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
		//printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
	}

    inline void insert( OBJECT o, const Vec3d& pos, double r ){
        Vec3d dpos; Vec3i ipos;
        pos2box( pos, ipos,dpos );
        insert( o, ipos, dpos, r );
    }

};

#endif

