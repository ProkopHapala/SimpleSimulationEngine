
#ifndef  SimplexGrid_h
#define  SimplexGrid_h

#include "HashMap.h"

typedef unsigned short  UHALF;
const UHALF MAP_OFFSET = 0x7FFF;

inline bool simplexIndex( double x, double y, int& ia, int& ib, double& da, double& db ){
	ia = (int)y * 0.86602540378;
	double b = x - 0.5d*y;
	ib = (int)b;
	da = y - ia;
    db = b - ib;
    return ( da + db ) > 1.0d;
};

template < class NODE, class TILE >
class SimplexField{
    public:
    NODE node;
    TILE lo;
    TILE hi;

    inline void setTile( TILE t, bool s ){ if( s ){ hi = t;    }{ lo  =  t;   } }
    inline TILE getTile(         bool s ){ if( s ){ return hi; }{ return lo ; } }

};

template < class OBJECT >
class SimplexGrid : public HashMap<OBJECT>{
    public:
    double        step, invStep;

    inline void simplexIndexBare( double x, double y, UHALF& ia, UHALF& ib ) const {
        double a = ( invStep *   y * 0.86602540378     ) + MAP_OFFSET;
        double b = ( invStep * ( x - 0.5d*y          ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
    }

    inline bool simplexIndex( double x, double y, UHALF& ia, UHALF& ib, double& da, double& db ) const {
        double a = ( invStep *   y * 0.86602540378    ) + MAP_OFFSET;
        double b = ( invStep * ( x - 0.5d*y         ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
        da = y - ia;
        db = b - ib;
        return ( da + db ) > 1.0d;
    }

    inline void nodePoint( UHALF ia, UHALF ib, double& x, double& y ) const {
        UHALF ia_ = ia            - MAP_OFFSET;
        y = step * 0.86602540378 *( ia_ );
        x = step *                ( ib + 0.5d*ia_ - MAP_OFFSET );
    };

    inline void tilePoint( UHALF ia, UHALF ib, bool s, double& x, double& y ) const {
        nodePoint( ia, ia, x, y );
        if( s ){ x+=step; y+=0.57735026919*step; }else{ x+=0.5d*step; y+=0.28867513459*step;  }
    };

    //inline ULONG  getBucketInt    ( UHALF  ix, UHALF iy )const{ return ( iy << 16 ) + ix;                     };
	//inline void   unfoldBucketInt ( ULONG bucket, UHALF& ix, UHALF& iy  )const{ ix = bucket&0xFFFF; iy = (bucket>>16)&0xFFFF; }

	inline ULONG  getBucketInt    (               UHALF  ia, UHALF ib   )const{ return ( ib << 16 ) + ia;                     };
	inline void   unfoldBucketInt ( ULONG bucket, UHALF& ia, UHALF& ib  )const{ ia = bucket&0xFFFF; ib = (bucket>>16)&0xFFFF; };
    inline ULONG  getBucket       ( double  x, double y )const{ UHALF ia,ib; simplexIndexBare( x, y, ia, ib ); printf( " getBucket %3.3f %3.3f %i %i\n", x,y,ia,ib ); return getBucketInt( ia, ib );    };

    inline void init( double step_, UINT power_ ){
        step    = step_;
		invStep = 1/step;
		HashMap<OBJECT>::init( power_ );
    }

	// TODO: many of this function can be STATIC
	inline int  findInt            ( OBJECT* p, UHALF ia, UHALF ib    )const{ return HashMap<OBJECT>::find            ( p, getBucketInt( ia,ib ) ); };
	inline int  find               ( OBJECT* p, double x, double y    )const{ return HashMap<OBJECT>::find            ( p, getBucket   (  x, y ) ); };
	inline int  insertNoTestInt    ( OBJECT* p, UHALF ia, UHALF ib    )     { return HashMap<OBJECT>::insertNoTest    ( p, getBucketInt( ia,ib ) ); };
	inline int  insertNoTest       ( OBJECT* p, double x, double y    )     { return HashMap<OBJECT>::insertNoTest    ( p, getBucket   (  x, y ) ); };
	inline int  insertIfNewInt     ( OBJECT* p, UHALF ia, UHALF ib    )     { return HashMap<OBJECT>::insertIfNew     ( p, getBucketInt( ia,ib ) ); };
	inline int  insertIfNew        ( OBJECT* p, double x, double y    )     { return HashMap<OBJECT>::insertIfNew     ( p, getBucket   (  x, y ) ); };
	inline bool tryRemoveInt       ( OBJECT* p, UHALF ia, UHALF ib    )     { return HashMap<OBJECT>::tryRemove       ( p, getBucketInt( ia,ib ) ); };
	inline bool tryRemove          ( OBJECT* p, double x, double y    )     { return HashMap<OBJECT>::tryRemove       ( p, getBucket   (  x, y ) ); };
	inline UINT getBucketIndexesInt( UHALF  ia, UHALF ib, UINT * outi )const{ return HashMap<OBJECT>::getBucketIndexes(    getBucketInt( ia,ib ), outi );  }
	inline UINT getBucketIndexes   ( double  x, double y, UINT * outi )const{ return HashMap<OBJECT>::getBucketIndexes(    getBucket   (  x, y ), outi );  }
	inline UINT getBucketObjectsInt( UHALF  ia, UHALF ib, OBJECT**out )const{ return HashMap<OBJECT>::getBucketObjects(    getBucketInt( ia,ib ), out  );  }
	inline UINT getBucketObjects   ( double  x, double y, OBJECT**out )const{ return HashMap<OBJECT>::getBucketObjects(    getBucket   (  x, y ), out  );  }

};

#endif

