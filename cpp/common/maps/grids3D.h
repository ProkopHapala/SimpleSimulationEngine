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

class GridRulerInterface{public:};

int constexpr NBIT3D   = (int)(64/3);
int constexpr OFFSET3D = 1<<(NBIT3D-1);
int constexpr XMASK3D = (1<<(NBIT3D))-1;
int constexpr YMASK3D = XMASK3D<<NBIT3D;
int constexpr ZMASK3D = YMASK3D<<NBIT3D;

class CubeGridRulerUnbound : public GridRulerInterface { public:
    // This should work for points in range ( -2^20 .. +2^20 ) == ( -1e+6 .. +1e+6 )
    double step;
    double invStep;

    inline void     ipcell  ( const Vec3d& pos, Vec3i& ipos ) const { ipos.x=(int)(pos.x*invStep); ipos.y=(int)(pos.y*invStep); ipos.z=(int)(pos.z*invStep); }
    inline uint64_t pos2long( const Vec3d& pos              ) const { Vec3i ip; ipcell(pos,ip); return ixyz2long(ip); };

    int overlap_BBox(const Vec3d& p0, const Vec3d& p1, uint64_t* icells ) const {
        Vec3i ip0; ipcell( p0, ip0 );
        Vec3i ip1; ipcell( p1, ip1 );
        ip1.x++; ip1.y++; ip1.z++;
        int n=0;
        for(int ix=ip0.x; ix<ip1.x; ix++){
            for(int iy=ip0.y; iy<ip1.y; iy++){
                for(int iz=ip0.z; iz<ip1.z; iz++){
                    icells[n]=ixyz2long({ix,iy,iz});
                    n++;
                }
            }
        }
        return n;
    }

    inline uint64_t ixyz2long( Vec3i      ip         ) const { return (ip.x+OFFSET3D) + ( ( (ip.y+OFFSET3D) + ((ip.z+OFFSET3D)<<NBIT3D) )<<NBIT3D );  }
    inline void     long2ixyz( uint64_t i, Vec3i& ip ) const { ip.x = i&XMASK3D; ip.y = (i&YMASK3D)>>NBIT3D; ip.z = (i&ZMASK3D)>>(NBIT3D*2); }

};

class CubeGridRuler : public GridRulerInterface { public:
    double step;
    double invStep;
    Vec3d  pos0;
    Vec3d  pmax;
    Vec3d  span;
    Vec3i  n;
    int    ntot,nxy;
    //int    nxy;
    inline void setn( Vec3i n_ ){ n = n_; nxy = n.x*n.y; ntot=nxy*n.z; }
    inline void setStep( double step_ ){ step=step_; invStep=1/step; };

    inline void setup( Vec3d pmin, Vec3d pmax_, double step ){
        setStep(step);
        pos0=pmin;
        pmax=pmax_;
        span=pmax-pos0;
        setn( { (int)(span.z*invStep+1), (int)(span.y*invStep+1), (int)(span.z*invStep+1) } );
    };

    inline void pos2box( const Vec3d& pos, Vec3i& ipos, Vec3d& dpos ) const {
        dpos.x = x2grid( pos.x-pos0.x, step, invStep, ipos.x );
        dpos.y = x2grid( pos.y-pos0.y, step, invStep, ipos.y );
        dpos.z = x2grid( pos.z-pos0.z, step, invStep, ipos.z );
        //printf( "(%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", pos.x, pos.y, pos.z, ipos.x,ipos.y,ipos.z);
    }

    inline Vec3i ipcell( const Vec3d& pos ) const { return (Vec3i){ (int)(pos.x-pos0.x)*invStep, (int)(pos.y-pos0.y)*invStep, (int)(pos.z-pos0.z)*invStep }; }
    inline int   icell ( const Vec3d& pos ) const { return ixyz2i(ipcell(pos)); }
    //inline int   icell ( const Vec3d& pos ) const { return ixyz2i({ (int)(pos.x-pos0.x)*invStep, (int)(pos.y-pos0.y)*invStep, (int)(pos.z-pos0.z)*invStep } ); }

    inline Vec3d box2pos( const Vec3i& ipos, const Vec3d& dpos ) const {
        return (Vec3d){
            step*ipos.x + pos0.x + dpos.x,
            step*ipos.y + pos0.y + dpos.y,
            step*ipos.z + pos0.z + dpos.z };
    }

    int overlap_Sphere( Vec3d pos, double r, int* icells ) const {
        int n=1;
        Vec3d dpos; Vec3i ipos;
        pos2box( pos, ipos,dpos );
        icells[0]=ixyz2i( ipos );   // insert
        //printf( "icells[0] %i (%i,%i,%i) (%f,%f,%f)\n", icells[0], ipos.x,ipos.y,ipos.z,   pos.x,pos.y,pos.z );
	    int    dix =0,diy =0,diz =0;
	    double dr2x=0,dr2y=0,dr2z=0;
	    double mr = step-r;
	    //printf( "(%f,%f,%f) %f %f \n", dpos.x,dpos.y,dpos.z, r, mr );
	    if     (  dpos.x < r  ){ icells[n]=ixyz2i( {ipos.x-1, ipos.y  , ipos.z  } ); n++;  dix=-1; dr2x = sq(     dpos.x); }
	    else if(  dpos.x > mr ){ icells[n]=ixyz2i( {ipos.x+1, ipos.y  , ipos.z  } ); n++;  dix=+1; dr2x = sq(step-dpos.x); }
	    if     (  dpos.y < r  ){ icells[n]=ixyz2i( {ipos.x  , ipos.y-1, ipos.z  } ); n++;  diy=-1; dr2y = sq(     dpos.y); }
	    else if(  dpos.y > mr ){ icells[n]=ixyz2i( {ipos.x  , ipos.y+1, ipos.z  } ); n++;  diy=+1; dr2y = sq(step-dpos.y); }
	    if     (  dpos.z < r  ){ icells[n]=ixyz2i( {ipos.x  , ipos.y  , ipos.z-1} ); n++;  diz=-1; dr2z = sq(     dpos.z); }
	    else if(  dpos.z > mr ){ icells[n]=ixyz2i( {ipos.x  , ipos.y  , ipos.z+1} ); n++;  diz=+1; dr2z = sq(step-dpos.z); }
	    double r2 = r*r;
	    if ( dr2x+dr2y      < r2 ){ icells[n]=ixyz2i( {ipos.x+dix, ipos.y+diy, ipos.z    } ); n++; }
	    if ( dr2x+dr2z      < r2 ){ icells[n]=ixyz2i( {ipos.x+dix, ipos.y    , ipos.z+diz} ); n++; }
	    if ( dr2y+dr2z      < r2 ){ icells[n]=ixyz2i( {ipos.x    , ipos.y+diy, ipos.z+diz} ); n++; }
	    if ( dr2x+dr2y+dr2z < r2 ){ icells[n]=ixyz2i( {ipos.x+dix, ipos.y+diy, ipos.z+diz} ); n++; }
        return n;
	    //if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
	    //printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
    }

    int overlap_BBox(const Vec3d& p0, const  Vec3d& p1, int* icells ) const {
        Vec3i ip0 = ipcell( p0 );
        Vec3i ip1 = ipcell( p1 );
        ip1.x++; ip1.y++; ip1.z++;
        int n=0;
        for(int ix=ip0.x; ix<ip1.x; ix++){
            for(int iy=ip0.y; iy<ip1.y; iy++){
                for(int iz=ip0.z; iz<ip1.z; iz++){
                    icells[n]=ixyz2i({ix,iy,iz});
                    n++;
                }
            }
        }
        return n;
    }

    int overlap_BBox(const Vec3d& p0, const  Vec3d& p1, uint64_t* icells ) const {
        Vec3i ip0 = ipcell( p0 );
        Vec3i ip1 = ipcell( p1 );
        ip1.x++; ip1.y++; ip1.z++;
        int n=0;
        for(int ix=ip0.x; ix<ip1.x; ix++){
            for(int iy=ip0.y; iy<ip1.y; iy++){
                for(int iz=ip0.z; iz<ip1.z; iz++){
                    icells[n]=ixyz2long({ix,iy,iz});
                    n++;
                }
            }
        }
        return n;
    }

    int overlap_Line(const Vec3d& p0, const Vec3d& p1, int* icells ){
    }

    int overlap_Triangle(const Vec3d& pa, const Vec3d& pb, const Vec3d& pc, int* icells ){
    }

    inline uint64_t ixyz2long( Vec3i      ip         ) const { return (ip.x+OFFSET3D) + ( ( (ip.y+OFFSET3D) + ((ip.z+OFFSET3D)<<NBIT3D) )<<NBIT3D );  }
    inline void     long2ixyz( uint64_t i, Vec3i& ip ) const { ip.x = i&XMASK3D; ip.y = (i&YMASK3D)>>NBIT3D; ip.z = (i&ZMASK3D)>>(NBIT3D*2); }

    inline int  ixyz2i( Vec3i ip         ) const { return ip.x + n.x*(ip.y + n.y*ip.z);          }
    inline void i2ixyz( int i, Vec3i& ip ) const { ip.z=i/nxy; i=i%nxy; ip.y=i/n.x; ip.x=i%n.x;  }

};

namespace Grids3D {

    // should we move this to Vec3.h ?
    inline void wrapIndex( Vec3i& ipos, const Vec3i& n ){
        ipos.x = wrap_index_fast( ipos.x, n.x );
        ipos.y = wrap_index_fast( ipos.y, n.y );
        ipos.z = wrap_index_fast( ipos.z, n.z );
    }


    //typedef KeyType  unit32_t;
    //typedef ValType  unit32_t;

    template<class OBJECT, void insert_func( int ix, int iy, int iz, int val) >
    inline void insert( OBJECT o, const CubeGridRuler& ruler, const Vec3d& pos, double r ){
        Vec3d dpos; Vec3i ipos;
        ruler.pos2box( pos, ipos,dpos );
        insert<insert_func>( o, ipos, dpos, r );
    }

    //   search 2x2x2 neighborhood of
    template<class OBJECT, void insert_func( int ix, int iy, int iz, int val) >
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
	    if ( dr2x+dr2y      < r2 ){ insert_func( ipos.x+dix, ipos.y+diy, ipos.z    , o ); }
	    if ( dr2x+dr2z      < r2 ){ insert_func( ipos.x+dix, ipos.y    , ipos.z+diz, o ); }
	    if ( dr2y+dr2z      < r2 ){ insert_func( ipos.x    , ipos.y+diy, ipos.z+diz, o ); }
	    if ( dr2x+dr2y+dr2z < r2 ){ insert_func( ipos.x+dix, ipos.y+diy, ipos.z+diz, o ); }

	    //if( (dix!=0)&&(diy!=0) ){ insert( o, ipos.x+dix, ipos.y+diy ); }
	    //printf( " %1.3f %1.3f  (%1.3f,%1.3f) (%i,%i) %1.3f \n", r, mr, dpos.x,dpos.y, dix, diy, dr2 );
    }



}

#endif
