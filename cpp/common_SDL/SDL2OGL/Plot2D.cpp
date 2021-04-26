
#include <SDL2/SDL_opengl.h>

#include <stdio.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include "Plot2D.h"  // THE HEADER


/////////////////////////////////
//  class        DataLine2D
/////////////////////////////////

void DataLine2D::update(){
    bounds.x0 = +1e+300;
    bounds.x1 = -1e+300;
    bounds.y0 = +1e+300;
    bounds.y1 = -1e+300;
    for(int i=0; i<n; i++){
        double x  = xs[i];
        double y  = ys[i];
        bounds.x0 = _min(bounds.x0,x);
        bounds.x1 = _max(bounds.x1,x);
        bounds.y0 = _min(bounds.y0,y);
        bounds.y1 = _max(bounds.y1,y);
        //printf( "DataLine2D::update[%i] %g(%g,%g) %g(%g,%g) \n", i, x,bounds.x0,bounds.x1, y,bounds.y0, bounds.y1 );
        //printf( "%i <%f..%f> <%f..%f> \n", i, bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
    }
    //printf( " DataLine2D::update <%f..%f> <%f..%f> \n", bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
};

void DataLine2D::draw(){
    if(xs==NULL) return;
    if(yfunc){ VecN::set(n,xs,yfunc,ys); }
    if(lineStyle =='-')Draw2D::plot( n, xs, ys );
    switch(pointStyle){
        case '.': Draw2D::plot_dots ( n, xs, ys ); break;
        case '+': Draw2D::plot_cross( n, xs, ys, pointSize    ); break;
        case 'x': Draw2D::plot_X    ( n, xs, ys, pointSize    ); break;
        case '3': Draw2D::plot_O    ( n, xs, ys, pointSize, 3 ); break;
        case '4': Draw2D::plot_O    ( n, xs, ys, pointSize, 4 ); break;
        case 'o': Draw2D::plot_O    ( n, xs, ys, pointSize, 8 ); break;
    }
}

int DataLine2D::render(){
    //printf( "DalaLine2D::render() \n" );
    if( glObj ) glDeleteLists(glObj,1);
    glObj = glGenLists(1);
    //printf( " 1 DataLine2D::render() %i \n", glObj  );
    glNewList(glObj, GL_COMPILE);
    draw();
    glEndList( );
    //printf( " 2 DataLine2D::render() %i \n", glObj  );
    return glObj;
};

void DataLine2D::view(){
    if(bView){
        Draw::setRGBA(clr);
        //printf( "DEBUG DataLine2D::view() %i \n", glObj  );
        glCallList( glObj );
        //printf("%i\n", glObj);
    }
};

DataLine2D::~DataLine2D(){
    if(xs &&(!bSharedX) ) delete [] xs;
    if(ys) delete [] ys;
    if( glObj ) glDeleteLists(glObj,1);
}

/////////////////////////////////
//  class      Plot2D
/////////////////////////////////

void Plot2D::update(){
    //printf( "Plot2D::update \n" );
    bounds.x0 = +1e+300;
    bounds.x1 = -1e+300;
    bounds.y0 = +1e+300;
    bounds.y1 = -1e+300;
    int i=0;
    for( DataLine2D * line : lines ){
        line->update();
        bounds.x0 = _min(bounds.x0, line->bounds.x0);
        bounds.x1 = _max(bounds.x1, line->bounds.x1);
        bounds.y0 = _min(bounds.y0, line->bounds.y0);
        bounds.y1 = _max(bounds.y1, line->bounds.y1);
        //printf( "Plot2D::update.l[%i]    <%f..%f> <%f..%f> \n", i, bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
        i++;
    }
    //printf( "Plot2D::update:  <%f..%f> <%f..%f> \n", bounds.x0, bounds.x1,  bounds.y0, bounds.y1 );
};

void Plot2D::autoAxes(double dx, double dy){
    //printf( "Plot2D::autoAxes \n" );
    int n0,n1;
    n0=(int)(bounds.x0/dx)-1;  axBounds.x0 = dx*n0;
    n1=(int)(bounds.x1/dx)+1;  axBounds.x1 = dx*n1; nXTicks=(n1-n0)+1;
    //printf("x(%g:%g) dx %g nx(%i:%i) nXTicks %i \n", bounds.x0, bounds.x1, dx,   n0, n1, nXTicks );
    n0=(int)(bounds.y0/dy)-1;  axBounds.y0 = dy*n0;
    n1=(int)(bounds.y1/dy)+1;  axBounds.y1 = dy*n1; nYTicks=(n1-n0)+1;
    //printf("y(%g:%g) dy %g ny(%i:%i) nYTicks %i \n", bounds.y0, bounds.y1, dy,   n0, n1, nYTicks );

    if( xTicks==NULL ) delete xTicks;
    if( yTicks==NULL ) delete yTicks;
    //double x0    = (axBounds.x0 - axPos.x);  x0 = 2*x0 - dx*(int)(x0/dx);
    //double y0    = (axBounds.y0 - axPos.y);  y0 = 2*y0 - dy*(int)(y0/dy);
    //printf("DEBUG 2.1.1\n");
    xTicks = new double[nXTicks]; VecN::arange( nXTicks, axBounds.x0, dx, xTicks );
    //printf("DEBUG 2.1.2 %i\n", nYTicks );
    yTicks = new double[nYTicks];
    //printf("DEBUG 2.1.3\n");
    VecN::arange( nYTicks, axBounds.y0, dy, yTicks );
    //printf("DEBUG 2.1.4\n");
}

void Plot2D::drawTexts(){
    char str[16];
    if(bTicks){
    Draw::setRGBA(clrTicksX);
    Draw2D::drawText( xlabel.c_str(), 0, {shift.x,shift.y-2*tickSz*scaling.y}, 0.0, fontTex, tickSz );
    for(int i=0; i<nXTicks; i++){
        if(logX) { sprintf(str,tickFormat,pow(10,xTicks[i])); }
        else     { sprintf(str,tickFormat,xTicks[i]); }
        Draw2D::drawText(str, 0, {shift.x+xTicks[i]*scaling.x,shift.y+axPos.y*scaling.y}, 90, fontTex, tickSz );
    }
    Draw::setRGBA(clrTicksY);
    Draw2D::drawText( ylabel.c_str(), 0, (Vec2d){shift.x,shift.y-tickSz*ylabel.length()*scaling.y}, 90.0, fontTex, tickSz );
    for(int i=0; i<nYTicks; i++){
        sprintf(str,tickFormat,yTicks[i]);
        if(logY) { sprintf(str,tickFormat,pow(10,yTicks[i])); }
        else     { sprintf(str,tickFormat,yTicks[i]); }
        Draw2D::drawText(str, 0, {shift.x+axPos.x*scaling.x,shift.x+yTicks[i]*scaling.y}, 0.0, fontTex, tickSz );
    }
    }
    for( DataLine2D* line : lines ){ // lines labels
        Draw::setRGBA(line->clr);
        //Draw2D::drawLine({5.0,i+0.5},{0.0,0.0});
        //Draw2D::drawText(line->label.c_str(), 0, (Vec2d)legend_pos+(Vec2d){0,tickSz*i}, 0., fontTex, tickSz );
        //Draw2D::drawText(line->label.c_str(), 0, {line->xs[0],line->ys[0]}, 0., fontTex, tickSz );
        Draw2D::drawText(line->label.c_str(), 0, {shift.x+line->xs[0]*scaling.x
                                                 ,shift.y+line->ys[0]*scaling.y}, 0.0, fontTex, tickSz );
        Draw2D::drawText(line->label.c_str(), 0, {shift.x+line->xs[line->n-1]*scaling.x
                                                 ,shift.y+line->ys[line->n-1]*scaling.y}, 0.0, fontTex, tickSz );
    }
}

void Plot2D::drawAxes(){
    //Draw2D::drawPointCross({0.0,0.0},100.0);
    //grid=false;
    if(bGrid){
        Draw::setRGBA(clrGrid);
        Draw2D::drawGrid( nXTicks, xTicks, axBounds.y0, axBounds.y1, true  );
        Draw2D::drawGrid( nYTicks, yTicks, axBounds.x0, axBounds.x1, false );
    }
    if(bAxes){
        Draw::setRGBA(clrTicksX); Draw2D::drawLine( {axPos   .x ,axBounds.y0}, {axPos   .x ,axBounds.y1} );
        Draw::setRGBA(clrTicksY); Draw2D::drawLine( {axBounds.x0,axPos   .y }, {axBounds.x1,axPos   .y } );
    }
    if(bTicks){
        Draw::setRGBA(clrTicksY); Draw2D::drawGrid( nXTicks, xTicks, axPos.x, axPos.x+(tickSz/(scaling.y)), true  );
        Draw::setRGBA(clrTicksX); Draw2D::drawGrid( nYTicks, yTicks, axPos.x, axPos.y+(tickSz/(scaling.x)), false );
    }
}


int Plot2D::render(){
    //void drawGrid( bounds.x0, bounds.y0, bounds.x0, bounds.x0, dx, dy );
    if( glObj ) glDeleteLists(glObj,1);
    glObj = glGenLists(1);
    glNewList(glObj, GL_COMPILE);
    if( (clrBg&0xFF000000) ){ Draw::setRGBA( clrBg ); Draw2D::drawRectangle_d( axBounds.a, axBounds.b, true ); }
    drawAxes();
    int i=0;
    glEndList( );
    //char str[256];
    for( DataLine2D* line : lines ){
        //printf( "render line[%i]\n", i );
        line->render();
        i++;
    }
    // TO DO :
    //if( tickCaption ){ }
    return glObj;
}

void Plot2D::drawHline ( double y ){ Draw2D::drawLine_d( {axBounds.x0, y}, {axBounds.x1, y} ); };
void Plot2D::drawVline ( double x ){ Draw2D::drawLine_d( {x, axBounds.y0}, {x, axBounds.y1} ); };
//void Plot2D::drawCursor( Vec2d p, double sz ){Draw2D::drawPointCross_d(p); };

void Plot2D::view(bool bAxes){
    glPushMatrix();
    glTranslatef(shift.x,shift.y,0.0);
    glScalef    (scaling.x,scaling.y,1.0);
    if (bAxes) glCallList( glObj );
    for( DataLine2D* line : lines ){ line->view(); }
    // printf( "%i \n", glObj);
    glPopMatrix();
    if(fontTex){ drawTexts(); }
}

/*
void Plot2D::xsharingLines(int nl, int np){
    double * xs = new double[np];
    for(int il=0; il<nl; il++){
        lines.push_back( new DataLine2D() );
        lines.back()->n=np;
        lines.back()->xs=xs;
        lines.back()->ys=new double[np];
    }
}
*/

void Plot2D::xsharingLines(int nl, int np, double xmin, double dx, uint32_t* colors, int ncol){
    //double * xs = new double[np];
    //VecN::arange (np,xmin,dx,xs);
    float dc = 1./nl;
    if(colors==0){ colors=Draw::colors_rainbow; ncol=Draw::ncolors; }
    uint32_t col=0;
    if(ncol>0){Draw::colorScale(0,ncol,colors);}else{ col=colors[0]; }
    lines.push_back( new DataLine2D(np,xmin,dx,col) );
    double* xs = lines.back()->xs;
    for(int il=1; il<nl; il++){
        if(ncol>0){Draw::colorScale(dc*il,ncol,colors);}else{ col=colors[il]; }
        lines.push_back( new DataLine2D(np,xs,col) );
    }
}

//void savetxt(const char* fname, int n=-1, int* ils=0 ){
void Plot2D::savetxt(const char* fname){
    FILE* fout = fopen(fname,"w");
    for(int ip=0; lines[0]->n; ip++){
        fprintf( fout, " %g ", lines[0]->xs[ip] );
        for(int il=0; il<lines.size(); il++){
            fprintf( fout, " %g ", lines[il]->ys[ip] );
        }
        printf("\n");
    }
    fclose(fout);
}

void Plot2D::init( ){
    //axBounds.set( {-10.0,-10.0},{10.0,10.0});
    bounds.set( {-10.0,-10.0},{10.0,10.0});
    axPos .set( 0.0,0.0  );
    //init        ( 1.0, 1.0 );
    autoAxes( 1.0, 1.0 );
}

void Plot2D::clear( bool bDeep ){
    if(bDeep){
        for(DataLine2D* line: lines ){
            delete line;
        }
    }
    lines.clear();
    init();
}

void Plot2D::erase( int i ){
    delete lines[i];
    lines.erase( lines.begin() + i );
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



