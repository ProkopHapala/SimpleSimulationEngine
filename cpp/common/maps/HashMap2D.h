
#ifndef  HashMap2D_h
#define  HashMap2D_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "HashMap.h"
//#include "Map2D.h"

typedef unsigned short  UHALF;
//typedef short  UHALF;

const UHALF MAP_OFFSET = 0x7FFF;






template <class OBJECT >
class HashMap2D : public HashMap<OBJECT> {
	public:
	double        step, invStep;

	//inline UHALF  getIx     ( double  x )          const{ return (UHALF)( ((short)(invStep * x)) + MAP_OFFSET ); }; // does not work for negative
	//inline UHALF  getIy     ( double  y )          const{ return (UHALF)( ((short)(invStep * y)) + MAP_OFFSET ); };
	inline UHALF  getIx     ( double  x )          const{ return (UHALF)( ( (invStep * x) + MAP_OFFSET ) ); };
	inline UHALF  getIy     ( double  y )          const{ return (UHALF)( ( (invStep * y) + MAP_OFFSET ) ); };
	inline double getX      ( UHALF  ix )          const{ return step * ( ix - MAP_OFFSET );            };
	inline double getY      ( UHALF  iy )          const{ return step * ( iy - MAP_OFFSET );            };
	inline ULONG  getBucket ( UHALF  ix, UHALF iy )const{ return ( iy << 16 ) + ix;                     };
	inline ULONG  getBucket ( double  x, double y )const{ return getBucket( getIx(x), getIy(y) );       };
	inline void   unfoldBucket( ULONG bucket, UHALF& ix, UHALF& iy  )const{ ix = bucket&0xFFFF; iy = (bucket>>16)&0xFFFF; }
	inline void   unfoldBucket( ULONG bucket, double& x, double& y  )const{ UHALF ix,iy;  unfoldBucket( bucket, ix, iy  ); x=getX(ix); y=getY(iy); }

    inline void init( double step_, UINT power_ ){
        step    = step_;
		invStep = 1/step;
		HashMap<OBJECT>::init( power_ );
    }

	// TODO: many of this function can be STATIC
	inline int  find            ( OBJECT* p, UHALF ix, UHALF iy    )const{ return HashMap<OBJECT>::find            ( p, getBucket( ix,iy ) ); };
	inline int  find            ( OBJECT* p, double x, double y    )const{ return HashMap<OBJECT>::find            ( p, getBucket(  x, y ) ); };
	inline int  insertNoTest    ( OBJECT* p, UHALF ix, UHALF iy    ){ return HashMap<OBJECT>::insertNoTest    ( p, getBucket( ix,iy ) ); };
	inline int  insertNoTest    ( OBJECT* p, double x, double y    ){ return HashMap<OBJECT>::insertNoTest    ( p, getBucket(  x, y ) ); };
	inline int  insertIfNew     ( OBJECT* p, UHALF ix, UHALF iy    ){ return HashMap<OBJECT>::insertIfNew     ( p, getBucket( ix,iy ) ); };
	inline int  insertIfNew     ( OBJECT* p, double x, double y    ){ return HashMap<OBJECT>::insertIfNew     ( p, getBucket(  x, y ) ); };
	inline bool tryRemove       ( OBJECT* p, UHALF ix, UHALF iy    ){ return HashMap<OBJECT>::tryRemove       ( p, getBucket( ix,iy ) ); };
	inline bool tryRemove       ( OBJECT* p, double x, double y    ){ return HashMap<OBJECT>::tryRemove       ( p, getBucket(  x, y ) ); };
	inline UINT getBucketIndexes( UHALF  ix, UHALF iy, UINT * outi )const{ return HashMap<OBJECT>::getBucketIndexes( getBucket( ix,iy ), outi );  }
	inline UINT getBucketIndexes( double  x, double y, UINT * outi )const{ return HashMap<OBJECT>::getBucketIndexes( getBucket(  x, y ), outi );  }
	inline UINT getBucketObjects( UHALF  ix, UHALF iy, OBJECT**out )const{ return HashMap<OBJECT>::getBucketObjects( getBucket( ix,iy ), out  );  }
	inline UINT getBucketObjects( double  x, double y, OBJECT**out )const{ return HashMap<OBJECT>::getBucketObjects( getBucket(  x, y ), out  );  }

    // this is most general and safest method for inserting 3D objects
/*
    inline int insert( OBJECT* p ){
        int n_inserts = 0;
        Rect2d bbox;
        p.boundingBox( &bbox );
        UHALF ix0 = getIx( bbox.x0 ); // TODO: bound check ?
        UHALF iy0 = getIy( bbox.y0 );
        UHALF ix1 = getIx( bbox.x1 );
        UHALF iy1 = getIy( bbox.y1 );
        for( UHALF iy=iy0; iy<=iy1; iy++ ){
            for( UHALF ix=ix0; ix<=ix1; ix++ ){
                n_inserts++;
                insert( p, ix, iy );
                //int i = getIndex( ix, iy );
                //tiles[i].push_back( p );
            }
        }
        return n_inserts;
    }
*/

	UINT getObjectsInRect( double x0, double y0, double x1, double y1, OBJECT**out ){
		UHALF ix0 = getIx( x0 );  UHALF iy0 = getIy( y0 );
		UHALF ix1 = getIx( x1 );  UHALF iy1 = getIy( y1 );
		UINT nfound = 0;
		for( UHALF iy = iy0; iy<=iy1; iy++ ){
			for( UHALF ix = ix0; ix<=ix1; ix++ ){
				UINT ni = getBucketObjects( ix, iy, out+nfound );
				nfound += ni;
			}
		}
		return nfound;
	}

};

#endif

