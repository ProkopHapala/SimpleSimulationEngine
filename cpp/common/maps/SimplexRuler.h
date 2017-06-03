
#ifndef  SimplexRuler_h
#define  SimplexRuler_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

// a = (1.0, 0.0)
// b = (0.5, 0.86602540378)

//const Vec2i kind2edge[4][2] = {  {{0,0},{0,1}}, {{0,0},{1,0}}, {{0,1},{1,0}}, {{0,0},{1,0}} };
const Vec2i kind2edge[4][2] = {  {{0,0},{1,0}}, {{0,0},{0,1}}, {{1,0},{0,1}}, {{0,0},{0,1}} };

inline Vec2d  toBaricentric  ( Vec2d p ){ return (Vec2d){ p.x-0.57735026919*p.y, p.y*1.15470053839 }; }
inline Vec2d  fromBaricentric( Vec2d p ){ return (Vec2d){ p.x+0.5*p.y, 0.86602540378*p.y           }; }
inline double trinagleInterp ( Vec2d dind, Vec3d hs ){ return hs.dot( {dind.a, dind.b, 1-dind.a-dind.b} ); }

inline Vec2d trinagleDeriv( Vec3d hs ){
    double dFda = hs.a - hs.c;
    double dFdb = hs.b - hs.c;
    double dFdx = dFda;
    double dFdy =  -0.57735026919*dFda + 1.15470053839d*dFdb;
    return (Vec2d){dFdx,dFdy};
    //return toBaricentric( {dFda,dFdb} );
    //return fromBaricentric( {dFda,dFdb} );
}

class SimplexRuler{
    public:
    double    MAP_OFFSET = 0.0;
    double    step=1.0d, invStep=1.0d;

    int        na,nb,ntot;

    Vec3d ray_p, ray_ip, ray_md;
    Vec2i ray_i;
    double ray_t;

    // --- inline functions

    inline Vec2i i2ip(int   i ) const { return {i%na,i/na};    } // https://stackoverflow.com/questions/7070346/c-best-way-to-get-integer-division-and-remainder
    inline int   ip2i(Vec2i ip) const { return (na*ip.b+ip.a); }

    inline Vec2i wrap_index( Vec2i ip)const{ return (Vec2i){ wrap_index_fast(ip.a,na),  wrap_index_fast(ip.b,nb) }; }
    inline Vec2i clam_index( Vec2i ip)const{ return (Vec2i){clamp_index_fast(ip.a,na), clamp_index_fast(ip.b,nb) }; }

    inline void simplexIndexBare( const Vec2d& p, Vec2i& ind ) const {
        double a = ( invStep * (  p.x - 0.57735026919*p.y   ) ) + MAP_OFFSET;
        double b = ( invStep *    p.y * 1.15470053839         ) + MAP_OFFSET;
    	ind.a = (int)a;
        ind.b = (int)b;
    }

    inline int simplexIndex( const Vec2d& p, Vec2i& ind, Vec2d& dind ) const {
        double a = ( invStep * (  p.x - 0.57735026919*p.y   ) ) + MAP_OFFSET;
        double b = ( invStep *    p.y * 1.15470053839       ) + MAP_OFFSET;
    	ind.a  = (int)a;
        ind.b  = (int)b;
        dind.a = a - ind.a;
        dind.b = b - ind.b;
        //printf( " %i \n", (dind.a+dind.b)>1.0d );
        //return (ip2i(ind)<<1) + (int)( (dind.a+dind.b)>1.0d );
        return (dind.a+dind.b)>1.0d;
    }

    inline void nodePoint( const Vec2i& ind, Vec2d& p ) const {
        double ib_ = ind.b-MAP_OFFSET;
        p.x = step * ( (ind.a-MAP_OFFSET) + 0.5*ib_ );
        p.y = step * ib_ * 0.86602540378;
    };

    inline void nodePoint( int i, Vec2d& p ) const {
        Vec2i ind = i2ip( i );
        double ib_ = ind.b-MAP_OFFSET;
        p.x = step * ( (ind.a-MAP_OFFSET) + 0.5*ib_ );
        p.y = step * ib_ * 0.86602540378;
    };

    inline void tilePoint( const Vec2i& ind, bool s, Vec2d& p ) const {
        nodePoint( ind, p );
        if( s ){ p.x+=step; p.y+=0.57735026919*step; }else{ p.x+=0.5d*step; p.y+=0.28867513459*step;  }
    };

    inline void tilePoint( int i, Vec2d& p ) const { tilePoint( i2ip(i>>1), i&1, p ); };

    inline double getValue( Vec2d p, double * hs )const {
        Vec2i ind; Vec2d dind;
        bool s = simplexIndex( p, ind, dind );
        int ia  = ip2i(wrap_index({ind.a+1,ind.b  }));
        int ib  = ip2i(wrap_index({ind.a  ,ind.b+1}));
        int ic;
        if( s ){ ic  = ip2i(wrap_index({ind.a+1,ind.b+1}));  dind={1-dind.y,1-dind.x}; }
        else   { ic  = ip2i(wrap_index(ind));                                          };
        //if( s ){ ic  = ip2i(wrap_index({ind.a+1,ind.b+1}));  dind={1-dind.y,1-dind.x};   nodePoint(wrap_index({ind.a+1,ind.b+1 }), p );  Draw2D::drawCircle_d( p, 0.25, 16, true ); }
        //else   { ic  = ip2i(wrap_index(ind));                                            nodePoint(wrap_index({ind.a,  ind.b   }), p );  Draw2D::drawCircle_d( p, 0.25, 16, true ); };
        //nodePoint(wrap_index({ind.a+1,ind.b   }), p );   Draw2D::drawCircle_d( p, 0.25, 16, true );
        //nodePoint(wrap_index({ind.a  ,ind.b+1 }), p );   Draw2D::drawCircle_d( p, 0.25, 16, true );
        //printf("(%i,%i) (%i,%i,%i)\n",ind.a,ind.b, ia, ib, ic );
        return trinagleInterp( dind, {hs[ia],hs[ib],hs[ic]} );
    }

    inline Vec2d getDeriv( Vec2d p, double * hs )const {
        Vec2i ind; Vec2d dind;
        bool s = simplexIndex( p, ind, dind );
        int ia  = ip2i(wrap_index({ind.a+1,ind.b  }));
        int ib  = ip2i(wrap_index({ind.a  ,ind.b+1}));
        int ic;
        if( s ){ ic  = ip2i(wrap_index({ind.a+1,ind.b+1}));  return trinagleDeriv( {hs[ib],hs[ia],hs[ic]} )*-invStep; }
        else   { ic  = ip2i(wrap_index(ind));                return trinagleDeriv( {hs[ia],hs[ib],hs[ic]} )* invStep; };
    }

    inline void setSize( int na_, int nb_ ){ na=na_; nb=nb_; ntot=na*nb; };
    inline void setStep( double step_     ){ step = step_; invStep = 1.0d/step; };

    int rayStart( Vec2d ray0, Vec2d hray ){
        ray_p.a  = hray.dot( { 0.0d        ,1.15470053839*invStep} );
        ray_p.b  = hray.dot( { 1.0d*invStep,0.57735026919*invStep} );
        ray_p.c  = hray.dot( {-1.0d*invStep,0.57735026919*invStep} );
        ray_md.a = ray0.dot( { 0.0d        ,1.15470053839*invStep} );
        ray_md.b = ray0.dot( { 1.0d*invStep,0.57735026919*invStep} );
        ray_md.c = ray0.dot( {-1.0d*invStep,0.57735026919*invStep} );
        if( ray_p.a < 0 ){ ray_p.a=-ray_p.a; ray_md.a = 1-ray_md.a; };
        if( ray_p.b < 0 ){ ray_p.b=-ray_p.b; ray_md.b = 1-ray_md.b; };
        if( ray_p.c < 0 ){ ray_p.c=-ray_p.c; ray_md.c = 1-ray_md.c; };
        int ia=(int)(ray_md.a + MAP_OFFSET); ray_md.a = 1-(ray_md.a - (ia - MAP_OFFSET) );
        int ib=(int)(ray_md.b + MAP_OFFSET); ray_md.b = 1-(ray_md.b - (ib - MAP_OFFSET) );
        int ic=(int)(ray_md.c + MAP_OFFSET); ray_md.c = 1-(ray_md.c - (ic - MAP_OFFSET) );
        ray_ip.set_inv( ray_p );
        ray_t = 0;
        simplexIndexBare( ray0, ray_i );
    }

    int rayStep(){
        Vec3d tm; tm.set_mul(ray_md, ray_ip);
        int edgeKind = -1;
        if( tm.a < tm.b ){
           if( tm.a < tm.c ){  // a min
                ray_t    += tm.a;
                ray_md.a   = 1; ray_md.b -= ray_p.b*tm.a; ray_md.c -= ray_p.c*tm.a;
                ray_i.b++;
                edgeKind = 0;
           }else{            // c min
                ray_t    += tm.c;
                ray_md.a -= ray_p.a*tm.c; ray_md.b -= ray_p.b*tm.c; ray_md.c = 1;
                ray_i.a++;
                edgeKind = 1;
           }
        }else{
           if( tm.b < tm.c ){  // b min
                ray_t    += tm.b;
                ray_md.a -= ray_p.a*tm.b; ray_md.b = 1; ray_md.c -= ray_p.c*tm.b;
                edgeKind = 2;
           }else{            // c min
                ray_t    += tm.c;
                ray_md.a -= ray_p.a*tm.c; ray_md.b -= ray_p.b*tm.c; ray_md.c = 1;
                ray_i.a++;
                edgeKind = 3;
           }
        }
        return edgeKind;
    }

    double rayView( Vec3d ray0, Vec3d hray, double * hs ){
        rayStart( {ray0.x,ray0.y}, {hray.x,hray.y} );
        double ot = 0;
        double oh = ray0.z;
        double og = getValue( {ray0.x,ray0.y}, hs );
        ray0.z   += og;
        for(int i=0; i<100; i++){
            rayStep();
            Vec3d  p = ray0 + hray*ray_t;
            double g = getValue( {p.x,p.y}, hs );
            double h = p.z-g;
            //double h = ray0.z+hray.z*ray_t - g;
            if( h<0 ){
                double  f = oh/(oh-h);
                double dt = ray_t-ot;
                ray_t     = ot + f*(ray_t-ot);
                return      og + f*(g-og);
            }
            oh = h;
            og = g;
            ot = ray_t;
        }
    }

};

#endif

