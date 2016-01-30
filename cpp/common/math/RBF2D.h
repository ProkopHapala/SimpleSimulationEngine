
#ifndef  RBF2D_h
#define  RBF2D_h

//#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

class RBF2D{
    public:
    Vec2d  pos;
    double rmax;
    double height;

/*
    int     ncliffs = 0;
    Line2d * cliffs = NULL;

    inline bool checkCliffs( const Vec2d& d )const{
        if( cliffs!=NULL ){
            for( int i=0; i<ncliffs; i++ ){
                double dist = cliffs[i].dir.dot( d );
                if( dist < 0 ) return true;
            }
        }
        return false;
    }
*/

    inline double getVal( const Vec2d& where ) const{
        Vec2d d; d.set_sub( where, pos );
        double r2    = d.norm2();
        double r2max = rmax*rmax;
        double mr2   = r2max - r2;
        if( mr2 > 0 ){
            //if ( checkCliffs( d ) ) return 0;
            mr2 /= r2max;
            double fr  = mr2*mr2;
            return height * fr;
        }else{
            return 0;
        }
    }

    inline bool getGrad( const Vec2d& where, Vec2d& grad ) const{
        Vec2d d; d.set_sub( where, pos );
        double r2 = d.norm2();
        double r2max = rmax*rmax;
        double mr2   = r2max - r2;
        if( mr2 > 0 ){
            //if ( checkCliffs( d ) ) return 0;
            mr2 /= r2max;
            grad.set_mul( d, mr2*height );
            return true;
        }else{
            return false;
        }
    }

    inline double eval( const Vec2d& where, Vec2d& grad ) const{
        Vec2d d; d.set_sub( where, pos );
        double r2 = d.norm2();
        double r2max = rmax*rmax;
        double mr2   = r2max - r2;
        if( mr2 > 0 ){
            //if ( checkCliffs( d ) ) return 0;
            mr2 /= r2max;
            double fr  = mr2*mr2;
            grad.set_mul( d, height * mr2 );
            return           height * fr   ;
        }else{
            return 0;
        }
    }

};

#endif

