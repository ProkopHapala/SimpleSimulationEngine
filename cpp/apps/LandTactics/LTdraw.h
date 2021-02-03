
#ifndef LTdraw_h
#define LTdraw_h

#include "LTWorld.h"
#include "Draw2D.h"

void cmapHeight(double g){

    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);

    glColor3f(g,g,g);

    //return g;
}

void plotSiteFittness( const LTsurrounding& sur, const LTUnit& u, const Vec2d& pos ){
    double E = sur.unitPosFittness( &u, pos );
    sprintf( strBuf, "%3.3f", E );
    Draw2D::drawText(strBuf, pos, {100.0,20.0}, default_font_texture, 0.5 );
}

void plotSurrounding( const LTsurrounding& sur, const Vec2d& pos ){
    glColor3f( 0.0,0.0,1.0 ); for( LTUnit* u : sur.coleagues     ){ Draw2D::drawLine_d( pos, u->pos ); }
    glColor3f( 1.0,0.0,0.0 ); for( LTUnit* u : sur.enemies       ){ Draw2D::drawLine_d( pos, u->pos ); }
    glColor3f( 0.5,1.0,0.0 ); for( LTLinearObject* l : sur.lobjs ){ Draw2D::drawLine_d( pos, (l->p1+l->p2)*0.5 ); }
    glColor3f( 0.0,1.0,0.5 ); for( LTStaticObject* o : sur.objs  ){ Draw2D::drawLine_d( pos, o->pos ); }
}

void drawPath( SimplexRuler& ruler, const Way& way ){
    glBegin(GL_LINE_STRIP);
    Vec2i oip=Vec2iZero;
    for( int i : way.path ){
        Vec2d p;
        Vec2i ip = ruler.i2ip( i );
        ruler.nodePoint( ip, p );
        Vec2i dip = ip-oip;
        if( (dip.x>1)||(dip.x<-1)||(dip.y>1)||(dip.y<-1) ){
            glEnd();
            glBegin(GL_LINE_STRIP);
        }
        oip=ip;
        glVertex3f(p.x, p.y, 0.0 );
    }
    glEnd();
}

void drawRiver( SimplexRuler& ruler, const River& river ){
    glBegin(GL_LINE_STRIP);
    Vec2i oip=Vec2iZero;
    for( int i : river.path ){
        Vec2d p;
        Vec2i ip = ruler.i2ip( i );
        ruler.nodePoint( ip, p );
        Vec2i dip = ip-oip;
        if( (dip.x>1)||(dip.x<-1)||(dip.y>1)||(dip.y<-1) ){
            glEnd();
            glBegin(GL_LINE_STRIP);
        }
        oip=ip;
        glVertex3f(p.x, p.y, 0.0 );
    }
    glEnd();
}


double polyLinIntercept( int n, double* xs, double* ys, double y0, double a, const Rect2d& camSpan, const Vec2d& sc ){
    double ox = xs[0];
    double oy = ys[0];
    //for( int i=0; i<n; i++ ){ Draw2D::drawPointCross_d( { camXmin + xs[i]*xsc, camYmin + ysc*ys[i] }, 5.0 );  };
    for(int i=1; i<n; i++){
        double x  = xs[i];
        double y  = ys[i];
        double h  = (y0 + a*x ) -  y;
        glColor3f(0.0f,0.0f,1.0f); Draw2D::drawPointCross_d( { camSpan.a.x + x*sc.x, camSpan.a.y + sc.y*y     }, 1.0 );
        glColor3f(1.0f,0.0f,0.0f); Draw2D::drawPointCross_d( { camSpan.a.x + x*sc.x, camSpan.a.y + sc.y*(y+h) }, 1.0 );
        glColor3f(1.0f,0.0f,1.0f); Draw2D::drawLine( { camSpan.a.x+x*sc.x, camSpan.a.y + sc.y*y     },{ camSpan.a.x+x*sc.x, camSpan.a.y + sc.y*(y+h) } );
        //printf( "... %i x=%f %f \n", i, x, oy );
        if(h<0){
            double oh = (y0 + a*ox) - oy;
            double f  = oh/(oh-h);
            //printf( "OOO %i x=%f dx=%f %f %f oh=%f \n", i, x, x-ox, a, oy, oh );
            x  = f*(x-ox) + ox;
            //printf( "<<< %i %f %f %f dx=%f   %f %f %f \n", i, y0, h, x, (x-ox),   oh, f, x );

            glColor3f(0.0f,1.0f,0.0f); Draw2D::drawPointCross_d( { camSpan.a.x+x*sc.x, camSpan.a.y + sc.y*(y0 + a*x ) }, 1.0 );
            return x;
        }
        ox=x; oy=y;
    }
    return -1e-8;
}

/*
void drawTerrain(){
}
*/

void drawMap( SimplexRuler& ruler, double* vals,  double vmin, double vmax, int ncol=Draw::ncolors, const uint32_t * colors=&Draw::colors_rainbow[0] ){
    /*
    int ix,iy; Vec2d p;
    ix= 0;iy= 0; ruler.nodePoint( {ix,iy}, p  ); printf( "[%i,%i]->(%g,%g)\n", ix,iy,p.x,p.y );
    ix= 0;iy=10; ruler.nodePoint( {ix,iy}, p  ); printf( "[%i,%i]->(%g,%g)\n", ix,iy,p.x,p.y );
    ix=10;iy= 0; ruler.nodePoint( {ix,iy}, p  ); printf( "[%i,%i]->(%g,%g)\n", ix,iy,p.x,p.y );
    ix=10;iy=10; ruler.nodePoint( {ix,iy}, p  ); printf( "[%i,%i]->(%g,%g)\n", ix,iy,p.x,p.y );
    exit(0);
    */
    double scc = 1./(vmax-vmin);
    for (int iy=0; iy<ruler.na-1; iy++){
        glBegin( GL_TRIANGLE_STRIP );
        //int ii = (i0.y+iy)*ruler.n + i0.x;
        for (int ix=0; ix<ruler.nb; ix++){
            Vec2d p,p_;
            float val,val_;

            //glColor3f(1.0,1.0,1.0);
            val = vals[ruler.ip2i({ix,iy})];
            Draw::colorScale( (val-vmin)*scc, ncol, colors );
            ruler.nodePoint( {ix,iy}, p  );
            glVertex3f( p.x, p.y, 0 );

            val_ = vals[ruler.ip2i({ix,iy+1})];
            Draw::colorScale( (val_-vmin)*scc, ncol, colors );
            ruler.nodePoint( {ix,iy+1}, p_  );
            glVertex3f( p_.x, p_.y, 0 );
            //if((ix%10==0)&&(iy%10==0))printf( "drawMap[%i,%i] %g(%g,%g) %g(%g,%g)\n", ix,iy,  val,p.x,p.y, val_, p_.x,p_.y );
        }
        glEnd();
    }
}


template<typename Func>
void drawMap( Vec2i ns, Vec2d p0, Vec2d dp, Func func, double vmin, double vmax, int ncol=Draw::ncolors, const uint32_t * colors=&Draw::colors_rainbow[0] ){
    double scc = 1./(vmax-vmin);
    for (int iy=0; iy<ns.y-1; iy++){
        glBegin( GL_TRIANGLE_STRIP );
        //int ii = (i0.y+iy)*ruler.n + i0.x;
        for (int ix=0; ix<ns.x; ix++){
            Vec2d p;
            float val;

            p = p0 + (Vec2d){dp.x*ix,dp.y*iy};
            val = func( p );
            Draw::colorScale( (val-vmin)*scc, ncol, colors );
            glVertex3f( p.x, p.y, 0 );

            p.y+=dp.y;
            val = func( p );
            Draw::colorScale( (val-vmin)*scc, ncol, colors );
            glVertex3f( p.x, p.y, 0 );
            //if((ix%10==0)&&(iy%10==0))printf( "drawMap[%i,%i] %g(%g,%g) %g(%g,%g)\n", ix,iy,  val,p.x,p.y, val_, p_.x,p_.y );
        }
        glEnd();
    }
}





void drawVisibilityIsolines( LTWorld& world, Vec2d ray0, int ndhs, int nDirs, double phiMin, double phiMax, double dhMin, double dhMax, double tmax, bool bFill=true ){
    //const int  ndhs = 5;
    //const int nDirs = 50;
    double dhs[ndhs];
    double ts [ndhs];
    Vec2d  hRays    [nDirs];
    double horizonts[nDirs*ndhs];
    // renerate slopes
    double ddhs = (dhMax-dhMin)/(ndhs-1);
    for(int i=0; i<ndhs ; i++ ){ dhs[i]=dhMin+ddhs*i; };
    // generate directions
    double dphi = (phiMax-phiMin)/(nDirs-1);
    for(int i=0; i<nDirs; i++ ){ double a = phiMin+dphi*i; hRays[i]=(Vec2d){cos(a),sin(a)}; };
    // raytrace terrain
    //double tmax = 200.0;
    for(int i=0; i<nDirs; i++ ){
        int nhit = world.ruler.rayHorizonts( ray0, hRays[i], world.ground, 10.0, ndhs, dhs, ts, tmax );
        for(int j=0; j<ndhs ; j++ ){
            int ij = j*nDirs + i;
            if( j<nhit ){ horizonts[ij]=ts[j];    }
            else        { horizonts[ij]=tmax+1.0; }
        }
    };
    // plot lines
    float cstep = 1.0f/(ndhs-1);
    if(bFill){
    	glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glDisable(GL_DEPTH_TEST);
        //int j = iline;
        //float c = cstep*j;
            for(int j=0; j<ndhs; j++ ){
            glBegin(GL_TRIANGLE_FAN);
                glVertex3f( ray0.x, ray0.y, 100.0);
                for(int i=0; i<nDirs; i++ ){
                    int ij = j*nDirs + i;
                    double t = horizonts[ij];
                    Vec2d p = ray0 + hRays[i]*t;
                    //if( t>tmax ){ glColor3f(0.0f,0.0f,0.0f); }else{ glColor3f(1.0f-c,0.0f,c); };
                    glColor4f(1.0,0.0,1.0,0.05);
                    glVertex3f( (float)p.x, (float)p.y, 100.0);
                }
            glEnd();
            }
    }else{ // Plot As Lines
        for(int j=0; j<ndhs; j++ ){
            float c = cstep*j;
            glBegin(GL_LINE_LOOP);
            //glBegin(GL_LINE_STRIP);
                for(int i=0; i<nDirs; i++ ){
                    int ij = j*nDirs + i;
                    double t = horizonts[ij];
                    Vec2d p = ray0 + hRays[i]*t;
                    if( t>tmax ){ glColor3f(0.0f,0.0f,0.0f); }else{ glColor3f(1.0f-c,0.0f,c); };
                    glVertex3f( (float)p.x, (float)p.y, 100.0);
                }
            glEnd();
        }
    }
}

#endif
