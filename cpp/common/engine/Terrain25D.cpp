
#include "Terrain25D.h" // THE HEADER
#include "spline_hermite.h"

// ====================================
// ============= Terrain25D ===========
// ====================================

double Terrain25D::eval(const Vec2d& pos, Vec2d& deriv ){
    constexpr double scx = 0.4;
    constexpr double scy = 0.15;
    double x  = scx*pos.x;
    double y  = scy*pos.y;
    double cx = cos(x);
    double sx = sin(x);
    double cy = cos(y);
    double sy = sin(y);
    deriv.set( scx*cx, scy*cy );
    //printf( "deriv (%3.3f,%3.3f)  (%3.3f,%3.3f) \n", deriv.x, deriv.y, x, y );
    return sx + sy;
};

double Terrain25D::ray( const Vec3d& hRay, const Vec3d& ray0, double tmax, Vec3d& normal ){
    constexpr double scx = 0.1;
    constexpr double scy = 0.2;
    double tstart,tend;
    boundingPlanes( hRay.z, ray0.z, tmax, tstart, tend );
    double tspan = (tend-tstart);
    Vec3d p; p.set_lincomb(1,ray0, tspan,hRay);
    double dt = fabs( (hRay.z+(scx+scy)*0.72) );
    Vec3d dRay; dRay.set_mul( hRay, dt );
    int n = tspan/dt;
    Vec2d deriv;
    double old_val = eval( {p.x, p.y}, deriv );
    for(int i=0; i<n; i++){
        p.add(dRay);
        double val = eval( {p.x, p.y}, deriv );
        if( val*old_val < 0 ){
            normal.set( deriv.x, deriv.y, -1 );
            return tstart+i*dt;
        }
    }
    return 1e+300;
};

// ============================================
// ============= Terrain25D_bicubic ===========
// ============================================

double Terrain25D_bicubic::eval( const Vec2d& pos, Vec2d& deriv ){
    Vec2d dipos; Vec2i ipos;
    ruler.pos2index( pos, dipos, ipos );
    ipos.x--; ipos.y--;

    if( (ipos.x<0)||(ipos.y<0)||(ipos.y>(ruler.n.x-3))||(ipos.y>(ruler.n.y-3)) ){
        //printf( "invalid index %i %i \n", ipos.x, ipos.y );
        deriv.x = 0.0; deriv.y = 0.0;  return 0.0;
    };

    double * h0 = heights + ruler.ip2i( ipos );
    double * h1 = h0 + ruler.n.x;
    double * h2 = h1 + ruler.n.x;
    double * h3 = h2 + ruler.n.x;
    //double val = Spline_Hermite::val2D<double>( dipos.x, dipos.y, h0, h1, h2, h3 );
    //printf( " pos (%g,%g) ipos (%i,%i) %i dpos (%g,%g) val %g \n", pos.x, pos.y, ipos.x, ipos.y, ruler.ip2i( ipos ), dipos.x, dipos.y, val );
    double val = Spline_Hermite::dval2D<double>( dipos.x, dipos.y, deriv.x, deriv.y, h0, h1, h2, h3 );
    deriv.mul( ruler.invStep );
    return val;
};

/*
int TerrainCubic::rayLine( Vec2d hdir, Vec2d p0, double hg0, double dr, double rmax, int ntg, double * tgs, Vec3d * poss ){
    double h0 = getVal( p0.x, p0.y )+hg0;
    double r  = 0.0d;
    for(int itg=0; itg<ntg; itg++){
        double dh = tgs[itg]*dr;
        double h  = h0 + dh*r;
        double hg = h  - getVal( p0.x+hdir.x*r, p0.y+hdir.y*r );
        while( r<rmax ){
            double r_  = r+dr;
            double h_  = h+dh;
            double hg_ = h_ - getVal( p0.x+hdir.x*r_, p0.y+hdir.y*r_ );
            if( hg_<0 ){
                r += dr*hg/(hg+hg_);
                poss[itg] = {p0.x+hdir.x*r,p0.y+hdir.y*r, h0 + dh*r };
                //printf( "%i %g %g (%g,%g,%g)\n", itg, dh, r, poss[itg].z, poss[itg].x, poss[itg].z );
                break;
            }else{
                r=r_;  h=h_;
            }
        }
    }
}
*/
