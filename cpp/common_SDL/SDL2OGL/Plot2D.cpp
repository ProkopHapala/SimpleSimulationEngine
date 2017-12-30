
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include "Plot2D.h"  // THE HEADER


/////////////////////////////////
//  class        DataLine2D
/////////////////////////////////

void DataLine2D::update(){
    bounds.x0 = +1e-300; bounds.x1 = -1e-300; bounds.y0 = +1e-300; bounds.y0 = -1e-300;
    for(int i=0; i<n; i++){
        double x  = xs[i];
        double y  = ys[i];
        bounds.x0 = _min(bounds.x0,x);  bounds.x1 = _max(bounds.x1,x);
        bounds.y0 = _min(bounds.y0,y);  bounds.y1 = _max(bounds.y1,y);
        //printf( "%i (%f,%f) (%f,%f) \n", i, bounds.x0, bounds.y0,  bounds.x1, bounds.y1 );
        //printf( "%i <%f..%f> <%f..%f> \n", i, bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
    }
    //printf( "  ---  <%f..%f> <%f..%f> \n", bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
};

void DataLine2D::draw(){
    if(xs==NULL) return;
    if(yfunc){ VecN::set(n,xs,yfunc,ys); }
    if(lineStyle =='-')Draw2D::plot( n, xs, ys );
    switch(pointStyle){
        case '+': Draw2D::plot_cross( n, xs, ys, pointSize    ); break;
        case 'x': Draw2D::plot_X    ( n, xs, ys, pointSize    ); break;
        case '3': Draw2D::plot_O    ( n, xs, ys, pointSize, 3 ); break;
        case '4': Draw2D::plot_O    ( n, xs, ys, pointSize, 4 ); break;
        case 'o': Draw2D::plot_O    ( n, xs, ys, pointSize, 8 ); break;
    }
}

void DataLine2D::render(){
    if( glObj ) glDeleteLists(glObj,1);
    glObj = glGenLists(1);
    glNewList(glObj, GL_COMPILE);
    draw();
    glEndList( );
};

void DataLine2D::view(){
    Draw::setRGBA(clr);
    glCallList( glObj );
    //printf("%i\n", glObj);
};

/////////////////////////////////
//  class      Plot2D
/////////////////////////////////

void Plot2D::update(){
    bounds.x0 = +1e-300; bounds.x1 = -1e-300; bounds.y0 = +1e-300; bounds.y0 = -1e-300;
    for( DataLine2D * line : lines ){
        line->update();
        bounds.x0 = _min(bounds.x0, line->bounds.x0); bounds.x1 = _max(bounds.x1, line->bounds.x1);
        bounds.y0 = _min(bounds.y0, line->bounds.y0); bounds.y1 = _max(bounds.y1, line->bounds.y1);
        //printf( "    <%f..%f> <%f..%f> \n", bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
    }
    //printf( "    <%f..%f> <%f..%f> \n", bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
};

void Plot2D::autoAxes(double dx, double dy){
    int n0,n1;
    n0=(int)(bounds.x0/dx)-1;  axBounds.x0 = dx*n0;
    n1=(int)(bounds.x1/dx)+1;  axBounds.x1 = dx*n1; nXTicks=(n1-n0)+1;
    n0=(int)(bounds.y0/dy)-1;  axBounds.y0 = dy*n0;
    n1=(int)(bounds.y1/dy)+1;  axBounds.y1 = dy*n1; nYTicks=(n1-n0)+1;
    if( xTicks==NULL ) delete xTicks;
    if( yTicks==NULL ) delete yTicks;
    //double x0    = (axBounds.x0 - axPos.x);  x0 = 2*x0 - dx*(int)(x0/dx);
    //double y0    = (axBounds.y0 - axPos.y);  y0 = 2*y0 - dy*(int)(y0/dy);
    xTicks = new double[nXTicks]; VecN::arange( nXTicks, axBounds.x0, dx, xTicks );
    yTicks = new double[nYTicks]; VecN::arange( nYTicks, axBounds.y0, dy, yTicks );

}

void Plot2D::drawAxes(){

    //Draw2D::drawPointCross({0.0,0.0},100.0);

    if(grid){
        Draw::setRGBA(clrGrid);
        Draw2D::drawGrid( nXTicks, xTicks, axBounds.y0, axBounds.y1, true  );
        Draw2D::drawGrid( nYTicks, yTicks, axBounds.x0, axBounds.x1, false );
    }else{
        Draw::setRGBA(clrTicksY); Draw2D::drawGrid( nXTicks, xTicks, axPos.x, axPos.x+tickSz,     true  );  Draw2D::drawLine( {axPos   .x ,axBounds.y0}, {axPos   .x ,axBounds.y1} );
        Draw::setRGBA(clrTicksX); Draw2D::drawGrid( nYTicks, yTicks, axPos.x, axPos.y+tickSz,     false );  Draw2D::drawLine( {axBounds.x0,axPos   .y }, {axBounds.x1,axPos   .y } );
    }

    if(fontTex){
        char str[16];
        Draw::setRGBA(clrTicksX);
        for(int i=0; i<nXTicks; i++){
            sprintf(str,tickFormat,xTicks[i]);
            Draw2D::drawText(str, 0, {xTicks[i],axPos.y}, 90, fontTex, tickSz );
        }
        Draw::setRGBA(clrTicksY);
        for(int i=0; i<nYTicks; i++){
            sprintf(str,tickFormat,yTicks[i]);
            Draw2D::drawText(str, 0, {axPos.x,yTicks[i]}, 0.0, fontTex, tickSz );
        }
    }

}

void Plot2D::render(){
    //void drawGrid( bounds.x0, bounds.y0, bounds.x0, bounds.x0, dx, dy );

    if( glObj ) glDeleteLists(glObj,1);
    glObj = glGenLists(1);
    glNewList(glObj, GL_COMPILE);
    if( (clrBg&0xFF000000) ){ Draw::setRGBA( clrBg ); Draw2D::drawRectangle_d( axBounds.a, axBounds.b, true ); }
    drawAxes();
    glEndList( );
    for( DataLine2D* line : lines ){ line->render(); }
    // TO DO :
    //if( tickCaption ){ }
};

void Plot2D::drawHline ( double y ){ Draw2D::drawLine_d( {axBounds.x0, y}, {axBounds.x1, y} ); };
void Plot2D::drawVline ( double x ){ Draw2D::drawLine_d( {x, axBounds.y0}, {x, axBounds.y1} ); };
//void Plot2D::drawCursor( Vec2d p, double sz ){Draw2D::drawPointCross_d(p); };

void Plot2D::view(){
    glCallList( glObj );
    for( DataLine2D* line : lines ){ line->view(); }
   // printf( "%i \n", glObj);
};

void Plot2D::init( ){
    axBounds.set( {-10.0,-10.0},{10.0,10.0});
    axPos   .set( 0.0,0.0  );
    //init        ( 1.0, 1.0 );
    autoAxes( 1.0, 1.0 );
}


/////////////////////////////////
//  class      QuePlot2D
/////////////////////////////////

void QuePlot2D::init( int n_, int nlines_ ){
    n=n_;
    nlines=nlines_;
    lColors = new uint32_t[nlines];
    data    = new double* [nlines];
    ts = new double[n];
    for(int iline=0; iline<nlines; iline++){
        data[iline]    = new double[n];
        lColors[iline] = hash_Wang( iline*5454+1545 );
    }
}


//int get_index(int )

void QuePlot2D::next(double t){
    int oip = ip;
    nsamp++;
    ip=wrap_index(ip+1);
    ts[ip] = t;
    for(int iline=0; iline<nlines; iline++){
        data[iline][ip] = data[iline][oip];
    };
};

void QuePlot2D::draw( bool xoff, bool yoff ){
    //int i0     = wrap_index(ip+1);
    //int nvalid = (nsamp<n)?nsamp:n;
    int ii0 = nsamp-n+1; if (ii0<0) ii0=0;
    int i0  = wrap_index( ii0 );
    //printf( "nvalid %i \n", nvalid );
    double t0  = ts[i0];
    for(int iline=0; iline<nlines; iline++){
        glBegin(GL_LINE_STRIP);
        Draw::setRGBA(lColors[iline]);
        double * dline = data[iline];
        double y0 = dline[i0];
        for(int ii=ii0; ii<nsamp; ii++){
            int i = wrap_index( ii );
            glVertex3f( ts[i]-t0, dline[i]-y0, 0.0);
        }
        glEnd();
    }
}

void QuePlot2D::drawTrj3D( Vec3i which ){
    int ii0 = nsamp-n+1; if (ii0<0) ii0=0;
    int i0  = wrap_index( ii0 );
    glBegin(GL_LINE_STRIP);
    double * xs = data[which.x];
    double * ys = data[which.y];
    double * zs = data[which.z];
    for(int ii=ii0; ii<nsamp; ii++){
        int i = wrap_index( ii );
        glVertex3f( xs[i], ys[i], zs[i] );
    }
    glEnd();
}

void QuePlot2D::drawTrj3DPoints( Vec3i which, double pointSize ){
    int ii0 = nsamp-n+1; if (ii0<0) ii0=0;
    int i0  = wrap_index( ii0 );

    double * xs = data[which.x];
    double * ys = data[which.y];
    double * zs = data[which.z];
    if( pointSize > 0 ){
        for(int ii=ii0; ii<nsamp; ii++){
            int i = wrap_index( ii );
            Draw3D::drawPointCross( {xs[i], ys[i], zs[i]}, pointSize );
        }
    }else{
        glBegin(GL_POINTS);
        for(int ii=ii0; ii<nsamp; ii++){
            int i = wrap_index( ii );
            glVertex3f( xs[i], ys[i], zs[i] );
        }
        glEnd();
    }
}



