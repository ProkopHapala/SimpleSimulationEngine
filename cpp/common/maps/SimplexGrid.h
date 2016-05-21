
#ifndef  SimplexGrid_h
#define  SimplexGrid_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"


#include "HashMap.h"

typedef unsigned short  UHALF;
const UHALF MAP_OFFSET = 0x7FFF;

inline bool simplexIndex( double x, double y, int& ia, int& ib, double& da, double& db ){
	ia = (int)y * 0.86602540378d;
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
    //Vec2d avec,bvec,cvec;
    //Vec2d ainv,binv,cinv;

    inline void simplexIndexBare( double x, double y, UHALF& ia, UHALF& ib ) const {
        double a = ( invStep *    y * 1.15470053839       ) + MAP_OFFSET;
        double b = ( invStep * (  x - 0.57735026919*y   ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
    }

    inline bool simplexIndex( double x, double y, UHALF& ia, UHALF& ib, double& da, double& db ) const {
        double a = ( invStep *    y * 1.15470053839       ) + MAP_OFFSET;
        double b = ( invStep * (  x - 0.57735026919*y   ) ) + MAP_OFFSET;
    	ia = (int)a;
        ib = (int)b;
        da = a - ia;
        db = b - ib;
        return ( ( da + db ) > 1.0d );
    }

    inline void nodePoint( UHALF ia, UHALF ib, double& x, double& y ) const {
        int ia_ = ((int)ia) - MAP_OFFSET;
        y = step * ia_ * 0.86602540378;
        x = step * ( ( ((int)ib) - MAP_OFFSET ) + 0.5*ia_ );
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
		invStep = 1.0d/step;
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

/*
	int raster_line( Vec2d dirHat, Vec2d pos0, Vec2d pos1, Vec2d * hits ){
	    double t0    = dirHat.dot( pos0 );
	    double t1    = dirHat.dot( pos1 );
	    //double dtmax = t1-t2;
        double pa,pb,pc, invPa,invPb,invPc;
        int    ia,ib,ic,i;
        pa = dirHat.dot({ 0,1               } );  invPa = 1/pa;   ia=(int)(t0*invPa + MAP_OFFSET);
        pb = dirHat.dot({ 0.86602540378,0.5d} );  invPb = 1/pb;   ib=(int)(t0*invPb + MAP_OFFSET);
        pc = dirHat.dot({-0.86602540378,0.5d} );  invPc = 1/pc;   ic=(int)(t0*invPc + MAP_OFFSET);
        printf( " t_1,2  %f %f   p_a,b,c %f %f %f  \n", t0, t1, pa, pb, pc );
        double t = t0;
        double da,db,dc, mta,mtb,mtc;
        da = t0*invPa;   ia=(int)(da + MAP_OFFSET);   da -= ( ia - MAP_OFFSET );   //tma = mda * invPa;
        db = t0*invPb;   ib=(int)(db + MAP_OFFSET);   db -= ( ib - MAP_OFFSET );   //tmb = mdb * invPb;
        dc = t0*invPc;   ic=(int)(dc + MAP_OFFSET);   dc -= ( ic - MAP_OFFSET );   //tmc = mdc * invPc;
        i=0;
        //exit(0);
        while( t<t1 ){
            double ta = da * invPa; double abs_ta = fabs( ta );
            double tb = db * invPb; double abs_tb = fabs( tb );
            double tc = dc * invPc; double abs_tc = fabs( tc );
            if( tma < tmb ){
               if( abs_ta < abs_tc ){  // a min
                    t    += abs_ta; ia++;
                    da    = 1; db -= pa*ta; dc -= pc*ta;
               }else{            // c min
                    t    += abs_tmc; ic++;
                    mda   -= pa*tmc; mdb -= pa*tmc; mdc = 1;
               }
            }else{
               if( abs_tmb < abs_tmc ){  // b min
                    t    += abs_tmb; ib++;
                    mda   = pa*tmb; mdb = 1; mdc -= pc*tmb;
                    ib++;
               }else{            // c min
                    t    += abs_tmc; ic++;
                    mda   -= pa*tmc; mdb -= pa*tmc; mdc = 1;
               }
            }
            hits[i].set_mul( dirHat, t );
            printf( "%i %f (%f,%f) \n", i, t, hits[i].x, hits[i].y );
            i++;
        }
        return i;
	}
*/




	int raster_line( Vec2d dirHat, Vec2d pos0, Vec2d pos1, Vec2d * hits, int * boundaries, UHALF * edges ){
	    double t0    = dirHat.dot( pos0 );
	    double t1    = dirHat.dot( pos1 );
	    double tspan = t1-t0;
        double pa,pb,pc, invPa,invPb,invPc;
        int    ia,ib,ic,i;
        printf( " %f %f \n", step, invStep );
        double mda,mdb,mdc, mta,mtb,mtc;
        pa  = dirHat.dot( { 0.0d        ,1.15470053839*invStep} );
        pb  = dirHat.dot( { 1.0d*invStep,0.57735026919*invStep} );
        pc  = dirHat.dot( {-1.0d*invStep,0.57735026919*invStep} );
        mda = pos0.dot  ( { 0.0d        ,1.15470053839*invStep} );
        mdb = pos0.dot  ( { 1.0d*invStep,0.57735026919*invStep} );
        mdc = pos0.dot  ( {-1.0d*invStep,0.57735026919*invStep} );
        if( pa < 0 ){ pa=-pa; mda = 1-mda; };
        if( pb < 0 ){ pa=-pb; mdb = 1-mdb; };
        if( pc < 0 ){ pc=-pc; mdc = 1-mdc; };
        ia=(int)(mda + MAP_OFFSET);   mda = 1-(mda - (ia - MAP_OFFSET) );
        ib=(int)(mdb + MAP_OFFSET);   mdb = 1-(mdb - (ib - MAP_OFFSET) );
        ic=(int)(mdc + MAP_OFFSET);   mdc = 1-(mdc - (ic - MAP_OFFSET) );
        invPa = 1/pa; invPb = 1/pb; invPc = 1/pc;
        //printf( " t_1,2  %f %f   p_a,b,c %f %f %f  \n", t0, t1, pa, pb, pc );
        printf( " pa invPa \n", pa, invPa );
        double t = 0;
        i=0;
        UHALF ia_,ib_;
        simplexIndexBare( pos0.x, pos0.y, ia_, ib_ );
        while( t<tspan ){
            double tma = mda * invPa;
            double tmb = mdb * invPb;
            double tmc = mdc * invPc;
            //t += tma; boundaries[i] = 0;  mda = 1;
            //t += tmb; boundaries[i] = 1;  mdb = 1;
            //t += tmc; boundaries[i] = 2;  mdc = 1;
            int ii = i<<2;
            if( tma < tmb ){
               if( tma < tmc ){  // a min
                    t    += tma;
                    mda   = 1; mdb -= pb*tma; mdc -= pc*tma;
                    boundaries[i] = 0; ia_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_; edges[ii+3] = ib_+1;
               }else{            // c min
                    t    += tmc;
                    mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                    boundaries[i] = 2; ib_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }
            }else{
               if( tmb < tmc ){  // b min
                    t    += tmb;
                    mda  -= pa*tmb; mdb = 1; mdc -= pc*tmb;
                    boundaries[i] = 1;
                    edges[ii  ] = ia_;   edges[ii+1] = ib_+1;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }else{            // c min
                    t    += tmc;
                    mda  -= pa*tmc; mdb -= pb*tmc; mdc = 1;
                    boundaries[i] = 2; ib_++;
                    edges[ii  ] = ia_; edges[ii+1] = ib_;
                    edges[ii+2] = ia_+1; edges[ii+3] = ib_;
               }
            }
            hits[i].set( pos0 );
            hits[i].add_mul( dirHat, t );
            //hits[i].set_mul( dirHat, t );
            //printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], tma, tmb, tmc,       t, hits[i].x, hits[i].y );
            printf( "%i %i  (%f,%f,%f)     %f (%f,%f) \n", i, boundaries[i], mda, mdb, mdc,       t, hits[i].x, hits[i].y );
            i++;
        }
        return i;
	}



};

#endif

