
#ifndef  TileBuffer_h
#define  TileBuffer_h

//template <,N>
//class BufferTile

#include "Vec2.h"

template<class OBJECT,int NX_,int NY_,int M_>
class TileBuffer{
    public:
    static constexpr int NX   = NX_;
	static constexpr int NY   = NY_;
    static constexpr int M    = M_;
	static constexpr int MNX  = M*NX;
	static constexpr int NXY  = NX*NY;
	static constexpr int NTOT = M*NX*NY;

	int    count;
	int    counts [NXY ];
	OBJECT buff   [NTOT];

	inline int  xy2i(        int  ix, int  iy ){ return ix + NX*iy; }
	inline void i2xy( int i, int& ix, int& iy ){ iy=i/NX; ix=i%NX;    } // see performance http://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder

	inline int  xym2i(        int  ix, int  iy, int  im ){ return im + M*ix + MNX*iy;         }
	inline void i2xym( int i, int& ix, int& iy, int& im ){ iy=i/MNX; i=i%MNX; ix=i/M; im=i%M; } // see performance http://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder

    inline int  im2i( int ixy, int im ){ return im + M*ixy; }

	void insert( OBJECT o, int ix, int iy ){
		int i = xy2i( ix, iy );
		count++;
		int im = counts[i];
		counts[i] = im+1;
		buff[ (i*M)+im ] = o;
	}

	int clear( int ix, int iy ){
		int i     = xy2i( ix, iy );
		int ni    = counts[i];
		counts[i] = 0;
		count    -= ni;
	}

	void clear( ){
		count = 0;
		for( int i=0; i<NXY; i++ ){
			counts[i] = 0;
		}
	}

	void insert( OBJECT o, const Vec2i& ipos, const Vec2d& dpos, double r ){
		insert( o, ipos.x, ipos.y  );
		int dix=0,diy=0;
		double dr2 = 0;
		double mr = 1-r;
		if     (  dpos.x < r  ){ insert( o, ipos.x-1, ipos.y );  dix=-1; dr2 += sq(  dpos.x); }
		else if(  dpos.x > mr ){ insert( o, ipos.x+1, ipos.y );  dix=+1; dr2 += sq(1-dpos.x); }
		if     (  dpos.y < r  ){ insert( o, ipos.x, ipos.y-1 );  diy=-1; dr2 += sq(  dpos.y); }
		else if(  dpos.y > mr ){ insert( o, ipos.x, ipos.y+1 );  diy=+1; dr2 += sq(1-dpos.y); }
		if ( dr2 < (r*r) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
		//if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
		//printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
	}

	OBJECT get( int ixy, int im ){
        return buff[ im2i( ixy, im ) ];
	}

};

#endif

