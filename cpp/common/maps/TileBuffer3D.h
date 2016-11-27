
#ifndef  TileBuffer3D_h
#define  TileBuffer3D_h

//template <,N>
//class BufferTile

#include "fastmath.h"
#include "Vec3.h"

template<class OBJECT,int NX_,int NY_,int NZ_,int M_>
class TileBuffer3D{
    public:
    static constexpr int NX   = NX_;
	static constexpr int NY   = NY_;
	static constexpr int NZ   = NZ_;
    static constexpr int M    = M_;
	static constexpr int MNX  = M*NX;
	static constexpr int NXY  = NX*NY;
	static constexpr int NXYZ = NX*NY*NZ;
	static constexpr int NTOT = M *NXYZ;

	Vec3d  pos0;
	double step,invStep;
	int    count;
	int    counts [NXYZ];
	OBJECT buff   [NTOT];

	inline int  xyz2i(        int  ix, int  iy, int  iz ){ return ix + NX*iy + NXY*iz; }
	inline void i2xyz( int i, int& ix, int& iy, int& iz ){ iz=i/NXY; i=i%NXY; iy=i/NX; ix=i%NX; } // see performance http://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder

    inline int  im2i( int ixyz, int im ){ return im + M*ixyz; };
	inline int  xyzm2i(        int  ix, int  iy, int  iz, int  im ){ return im2i( xyz2i(ix,iy,iz), im);  }
	inline void i2xyzm( int i, int& ix, int& iy, int& iz, int& im ){ im=i%M; i2xyz( i/M, ix, iy, iz);    }


	void insert( OBJECT o, int ix, int iy,  int iz ){
		int i = xyz2i( ix, iy, iz );
		count++;
		int im = counts[i];
		counts[i] = im+1;
		buff[ (i*M)+im ] = o;
	}

	int clear( int ix, int iy, int iz ){
		int i     = xyz2i( ix, iy, iz );
		int ni    = counts[i];
		counts[i] = 0;
		count    -= ni;
	}

	void clear( ){
		count = 0;
		for( int i=0; i<NXYZ; i++ ){
			counts[i] = 0;
		}
	}

	void insert( OBJECT o, const Vec3i& ipos, const Vec3d& dpos, double r ){
		insert( o, ipos.x, ipos.y, ipos.z );
		int dix=0,diy=0,diz=0;
		double dr2x=0,dr2y=0,dr2z=0;
		double mr = 1-r;
		if     (  dpos.x < r  ){ insert( o, ipos.x-1, ipos.y  , ipos.z  );  dix=-1; dr2x = sq(  dpos.x); }
		else if(  dpos.x > mr ){ insert( o, ipos.x+1, ipos.y  , ipos.z  );  dix=+1; dr2x = sq(1-dpos.x); }
		if     (  dpos.y < r  ){ insert( o, ipos.x  , ipos.y-1, ipos.z  );  diy=-1; dr2y = sq(  dpos.y); }
		else if(  dpos.y > mr ){ insert( o, ipos.x  , ipos.y+1, ipos.z  );  diy=+1; dr2y = sq(1-dpos.y); }
		if     (  dpos.z < r  ){ insert( o, ipos.x  , ipos.y  , ipos.z-1);  diz=-1; dr2z = sq(  dpos.z); }
		else if(  dpos.z > mr ){ insert( o, ipos.x  , ipos.y  , ipos.z+1);  diz=+1; dr2z = sq(1-dpos.z); }

		double r2 = r*r;
		if ( dr2x+dr2y      < r2 ){ insert( o, ipos.x+dix, ipos.y+diy, ipos.z     ); }
		if ( dr2x+dr2z      < r2 ){ insert( o, ipos.x+dix, ipos.y    , ipos.z+diz ); }
		if ( dr2y+dr2z      < r2 ){ insert( o, ipos.x    , ipos.y+diy, ipos.z+diz ); }
		if ( dr2x+dr2y+dr2z < r2 ){ insert( o, ipos.x+dix, ipos.y+diy, ipos.z+diz ); }

		//if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
		//printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
	}

    inline void pos2box( const Vec3d& pos, Vec3i& ipos, Vec3d& dpos ){
        dpos.x = x2grid( pos.x-pos0.x, step, invStep, ipos.x );
        dpos.y = x2grid( pos.y-pos0.y, step, invStep, ipos.y );
        dpos.z = x2grid( pos.z-pos0.z, step, invStep, ipos.z );
        //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", pos.x, pos.y, pos.z, ipos.x,ipos.y,ipos.z);
    }

    inline bool validIndex( const Vec3i& ipos ){
        return (ipos.x>=0)&&(ipos.y>=0)&&(ipos.z>=0)&&(ipos.x<NX)&&(ipos.y<NY)&&(ipos.z<NZ);
    }

    inline void insert( OBJECT o, const Vec3d& pos, double r ){
        Vec3d dpos; Vec3i ipos;
        pos2box( pos, ipos,dpos );
        insert( o, ipos, dpos, r );
    }

	inline OBJECT get( int ixyz, int im ){
        return buff[ im2i( ixyz, im ) ];
	}

};

#endif

