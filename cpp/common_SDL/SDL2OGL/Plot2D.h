
#ifndef  Plot2D_h
#define  Plot2D_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "VecN.h"

const int nXTicks_def = 11;
const int nYTicks_def = 11;
const int dXTicks_def = 1.0;
const int dYTicks_def = 1.0;

class DataLine2D{ public:
    // data
    int      n     =0;
    double * xs    =NULL;
    double * ys    =NULL;
    Func1d   yfunc =NULL;
    //update_data    =false;
    // draw properties
    Rect2d   bounds;
    int      glObj      =0;
    char     lineStyle  ='-';
    //char     pointStyle ='+';
    char     pointStyle =' ';
    float    pointSize  =0.1;
    uint32_t clr        =0xFFFF00FF;


    void update();
    void render();
    void draw();
    void view();

    inline void set     ( int n_, double * xs_, double * ys_){ n=n_; xs=xs_; ys=ys_; };
    inline void allocate( int n_ ){ n=n_; xs=new double[n]; ys=new double[n]; };

    inline void linspan(double xmin, double xmax){ VecN::linspan(n,xmin,xmax,xs); };

    inline DataLine2D();
    inline DataLine2D(int n_){ allocate(n_); }
};

class Plot2D{ public:
    // data
    std::vector<DataLine2D*> lines;
    // properties
    Rect2d   bounds;
    int      nXTicks=0,nYTicks=0;
    double * xTicks=NULL;
    double * yTicks=NULL;
    int      glObj=0;

    Vec2d  axPos;
    Rect2d axBounds;
    float tickSz = 0.1;

    bool     grid       =true;
    bool     tickCaption=false;
    uint32_t clrBg      = 0x00f0f0f0;
    uint32_t clrGrid    = 0xFFE0E0E0;
    uint32_t clrTicksX  = 0xFF000000;
    uint32_t clrTicksY  = 0xFF000000;
    char   * tickFormat = "%2.2f\0";
    int      fontTex=0;

    void update();
    void drawAxes();
    void render();
    void view  ();
    void init  ();
    //void init( double dx, double dy );
    void autoAxes(double dx, double dy);

    void drawHline ( double y );
    void drawVline ( double x );
    //void drawCursor( Vec2d p, double sz );
};

class QuePlot2D{ public:
    int n;
    int nlines;
    uint32_t * lColors;
    double   * ts;
    double   ** data;

    int nsamp=0;
    int ip=0;

    void init( int n_, int nlines_ );
    void next(double t);
    void draw( bool xoff, bool yoff );
    void drawTrj3D( Vec3i which );
    void drawTrj3DPoints( Vec3i which, double pointSize );

    inline int wrap_index(int i){ return  i>=n? i%n : i; };
    inline void set_back(int iline, double y ){ data[iline][ip] = y; };
};


#endif

