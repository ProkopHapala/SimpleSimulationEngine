
#ifndef LTdraw_h
#define LTdraw_h


void cmapHeight(double g){

    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);

    glColor3f(g,g,g);

    //return g;
}

/*
void drawStaticObject( LTStaticObject& o ){
    glBegin(GL_LINE_LOOP);
    glVertex3f( );
    glVertex3f( );
    glVertex3f( );
    glVertex3f( );
    glEnd();
    Draw2D::drawShape( o.pos, o.rot, o.type->glo );
}
*/

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


void drawTerrain(){


}



void drawVisibilityIsolines( LTWorld& world, Vec2d ray0, int ndhs, int nDirs, double phiMin, double phiMax, double dhMin, double dhMax, double tmax ){
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
    for(int j=0; j<ndhs; j++ ){
        //glBegin(GL_LINE_LOOP);
        float c = cstep*j;
        glBegin(GL_LINE_STRIP);
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

#endif
