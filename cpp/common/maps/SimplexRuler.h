
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
    double    step=1.0d, invStep=1.0d;
    int        na,nb,ntot;

    double    MAP_OFFSET    = 0.0;
    double    ROUND_OFFSET  = 1000.0; // WARRNING : if ROUND_OFFSET is insufficient rounding will fail !!!

    Vec3d  ray_d, ray_id, ray_p;
    Vec3i  ray_i,ray_si;
    double ray_t;

    int maxRayIter = 100;

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
        // dlat/dt
        ray_d.a = hray.dot( { 0.0d        ,1.15470053839*invStep} );
        ray_d.b = hray.dot( { 1.0d*invStep,0.57735026919*invStep} );
        ray_d.c = hray.dot( {-1.0d*invStep,0.57735026919*invStep} );
        // pos in grid coords
        ray_p.a = ray0.dot( { 0.0d        ,1.15470053839*invStep} );
        ray_p.b = ray0.dot( { 1.0d*invStep,0.57735026919*invStep} );
        ray_p.c = ray0.dot( {-1.0d*invStep,0.57735026919*invStep} );
        //printf( " d (%g,%g,%g) p (%g,%g,%g) \n", ray_d.a, ray_d.b, ray_d.c, ray_p.a, ray_p.b, ray_p.c );
        if( ray_d.a < 0 ){ ray_d.a=-ray_d.a; ray_p.a = 1-ray_p.a; ray_si.a=-1; }else{ ray_si.a=1; };
        if( ray_d.b < 0 ){ ray_d.b=-ray_d.b; ray_p.b = 1-ray_p.b; ray_si.b=-1; }else{ ray_si.b=1; };
        if( ray_d.c < 0 ){ ray_d.c=-ray_d.c; ray_p.c = 1-ray_p.c; ray_si.c=-1; }else{ ray_si.c=1; };
        //printf( " d (%g,%g,%g) p (%g,%g,%g) \n", ray_d.x, ray_d.y, ray_d.z, ray_p.x, ray_p.y, ray_p.z );
        // WARRNING : if ROUND_OFFSET is sufficient rounding will fail !!!
        ray_i.a=(int)(ray_p.a + ROUND_OFFSET); ray_p.a = 1-(ray_p.a - (ray_i.a - ROUND_OFFSET) );
        ray_i.b=(int)(ray_p.b + ROUND_OFFSET); ray_p.b = 1-(ray_p.b - (ray_i.b - ROUND_OFFSET) );
        ray_i.c=(int)(ray_p.c + ROUND_OFFSET); ray_p.c = 1-(ray_p.c - (ray_i.c - ROUND_OFFSET) );
        ray_id.set_inv( ray_d );
        ray_t = 0;
        //simplexIndexBare( ray0, ray_i );
        //printf( "p (%g,%g,%g) ip (%g,%g,%g)\n", ray_d.a, ray_d.b, ray_d.c,   ray_id.a, ray_id.b, ray_id.c );
    }

    int rayStep(){
        Vec3d tm; tm.set_mul(ray_p, ray_id);
        //printf( "tm (%g,%g,%g)  (%g,%g,%g) (%g,%g,%g)", tm.x, tm.y, tm.z, ray_md.x,    ray_md.y, ray_md.z,   ray_ip.x, ray_ip.y, ray_ip.z );
        //printf( "tm (%g,%g,%g) md (%g,%g,%g) \n", tm.a, tm.b, tm.c, ray_p.a, ray_p.b, ray_p.c,   ray_id.a, ray_id.b, ray_id.c );
        int edgeKind = -1;
        if( tm.a < tm.b ){
           if( tm.a < tm.c ){  // printf("a\n");
                ray_t   += tm.a;
                ray_p.a  = 1;
                ray_p.b -= tm.a*ray_d.b;
                ray_p.c -= tm.a*ray_d.c;
                ray_i.a +=ray_si.a;
                edgeKind = 0;
           }else{              // printf("c\n");
                ray_t   += tm.c;
                ray_p.a -= tm.c*ray_d.a;
                ray_p.b -= tm.c*ray_d.b;
                ray_p.c  = 1;
                ray_i.c +=ray_si.c;
                edgeKind = 1;
           }
        }else{
           if( tm.b < tm.c ){  // printf("b\n");
                ray_t   += tm.b;
                ray_p.a -= tm.b*ray_d.a;
                ray_p.b  = 1;
                ray_p.c -= tm.b*ray_d.c;
                ray_i.b +=ray_si.b;
                edgeKind = 2;
           }else{            // printf("c\n");
                ray_t   += tm.c;
                ray_p.a -= tm.c*ray_d.a;
                ray_p.b -= tm.c*ray_d.b;
                ray_p.c  = 1;
                ray_i.c +=ray_si.c;
                edgeKind = 1;
           }
        }
        return edgeKind;
    }

    double rayView( Vec2d ray0, Vec2d hray, double h0, double dh, double * hs, double tmax ){
        rayStart( ray0, hray );
        double ot = 0;
        double oh = h0;
        double og = getValue( ray0, hs );
        double z0 = og + h0;
        for(int i=0; i<maxRayIter; i++){
            rayStep();
            glColor3f(1.0f,0.0f,0.0f); Draw2D::drawPointCross_d( ray0 + hray*ray_t, 2.0 );
            Vec2d  p = ray0 + hray*ray_t;
            double z = z0   + dh  *ray_t;
            double g = getValue( p, hs );
            double h = z-g;
            Draw2D::drawLine_d( p, p+((Vec2d){-hray.y,hray.x})*h );
            //double h = ray0.z+hray.z*ray_t - g;
            if( h<0 ){
                double  f = oh/(oh-h);
                ray_t     = ot + f*(ray_t-ot);
                if(true){
                    // since the surface is piecewise-linear the solution must be exact !
                    double HH = (z0+dh*ray_t) - getValue( ray0+hray*ray_t, hs );
                    printf("HHH=%g ,   %g %g %g %g \n", HH, (oh-h), oh, h, f );
                    return HH;
                }
                return      og + f*(g-og);
            }else if( ray_t > tmax ){
                double  f = (tmax-ot)/(ray_t-ot);
                ray_t = tmax;
                return      og + f*(g-og);
            }
            oh = h;
            og = g;
            ot = ray_t;
        }
    }

    int rayList( Vec2d ray0, Vec2d hray, double h0, int ntg, double * tgs, double * Ttgs, double * hs, double tmax ){
        rayStart( ray0, hray );
        double ot = 0;
        //double oh = h0;
        //double og = getValue( ray0, hs );
        double z0 = getValue( ray0, hs ) + h0;
        for( int itg=0; itg<ntg; itg++){
            double dh = tgs[itg];
            double oh = (z0+dh*ot) - getValue( ray0+hray*ot, hs );
            if(oh<0){ printf("__error__ oh=%g itg \n", oh, itg ); }
            bool fly=true;
            while( fly ){
                // FIXME : new ray with "dh[n+1]" may still hit the same tile as dh[n] => we should not do rayStep() each time
                rayStep();
                Vec2d  p = ray0 + hray*ray_t;
                double z = z0   + dh  *ray_t;
                double g = getValue( p, hs );
                double h = z-g;
                if(oh<0){ printf("error oh=%g \n", oh ); }
                glColor3f(0.0f,1.0f,0.0f); Draw2D::drawLine_d( p, p+((Vec2d){-hray.y,hray.x})*h );
                if( h<0 ){
                    double  f = oh/(oh-h);
                    ray_t     = ot + f*(ray_t-ot);
                    if(f<0){printf("error f=%g (%g,%g) %i %g \n", f, oh, h,  itg, dh ); return itg; }
                    if(true){
                        double HH = (z0+dh*ray_t) - getValue( ray0+hray*ray_t, hs );
                        printf("HH=%g ,   %g %g %g \n", HH, (oh-h), oh, h );
                        //return itg;
                    }
                    fly = false;
                }else if( ray_t > tmax ){
                    ray_t     = tmax;
                    ot        = ray_t;
                    Ttgs[itg] = ray_t;
                    //fly = false;
                    return itg;
                }
                oh = h;
                ot = ray_t;
            }
            glColor3f(1.0f,0.0f,0.0f); Draw2D::drawPointCross_d( ray0 + hray*ray_t, 10.0 );
            Ttgs[itg] = ray_t;
        }
        return ntg;
    }

};

#endif

