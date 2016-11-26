
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"

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
    }
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
        bounds.x0 = _min(bounds.x0, line->bounds.x0); bounds.x1 = _min(bounds.x1, line->bounds.x1);
        bounds.y0 = _min(bounds.y0, line->bounds.y0); bounds.y1 = _min(bounds.y1, line->bounds.y1);
    }
};

void Plot2D::drawAxes(){

    //Draw2D::drawPointCross({0.0,0.0},100.0);

    if(grid){
        Draw::setRGBA(clrGrid);   Draw2D::drawGrid( nXTicks, xTicks, axBounds.y0, axBounds.y1, true  );
        Draw::setRGBA(clrGrid);   Draw2D::drawGrid( nYTicks, yTicks, axBounds.x0, axBounds.x1, false );
    }else{
        Draw::setRGBA(clrTicksY); Draw2D::drawGrid( nXTicks, xTicks, axPos.x, axPos.x+tickSz,     true  );  Draw2D::drawLine( {axPos   .x ,axBounds.y0}, {axPos   .x ,axBounds.y1} );
        Draw::setRGBA(clrTicksX); Draw2D::drawGrid( nYTicks, yTicks, axPos.x, axPos.y+tickSz,     false );  Draw2D::drawLine( {axBounds.x0,axPos   .y }, {axBounds.x1,axPos   .y } );
    }

}

void Plot2D::render(){
    //void drawGrid( bounds.x0, bounds.y0, bounds.x0, bounds.x0, dx, dy );

    if( glObj ) glDeleteLists(glObj,1);
    glObj = glGenLists(1);
    glNewList(glObj, GL_COMPILE);
    drawAxes();
    glEndList( );

    for( DataLine2D* line : lines ){ line->render(); }
    // TO DO :
    //if( tickCaption ){ }
};

void Plot2D::view(){
    for( DataLine2D* line : lines ){ line->view(); }
    glCallList( glObj );
   // printf( "%i \n", glObj);
};

void Plot2D::init( ){
    axBounds.set( {-10.0,-10.0},{10.0,10.0});
    axPos   .set( 0.0,0.0  );
    init        ( 1.0, 1.0 );
}

void Plot2D::init( double dx, double dy ){
    double xspan = axBounds.x1-axBounds.x0;
    double yspan = axBounds.y1-axBounds.y0;
    double x0    = (axBounds.x0 - axPos.x);  x0 = 2*x0 - dx*(int)(x0/dx);
    double y0    = (axBounds.y0 - axPos.y);  y0 = 2*y0 - dy*(int)(y0/dy);
    nXTicks=(int)((xspan)/dx)+1;
    nYTicks=(int)((yspan)/dy)+1;
    xTicks = new double[nXTicks]; VecN::arange( nXTicks, x0, dx, xTicks );
    yTicks = new double[nYTicks]; VecN::arange( nYTicks, y0, dx, yTicks );
    printf( "%i %i  %f %f\n", nXTicks, nYTicks, x0, y0 );
    printf( "%f %f  %f %f  %f %f\n", axBounds.x0, axBounds.y0, axBounds.x1, axBounds.y1, axPos.x, axPos.y );
    //VecN::print_vector( nXTicks, xTicks );
    //VecN::print_vector( nYTicks, yTicks );
    //exit(0);
};
