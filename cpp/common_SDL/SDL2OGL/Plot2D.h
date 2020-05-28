
#ifndef  Plot2D_h
#define  Plot2D_h

#include <vector>
#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "VecN.h"

const int nXTicks_def = 11;
const int nYTicks_def = 11;
const int dXTicks_def = 1.0;
const int dYTicks_def = 1.0;

// =======================
// ====  DataLine2D
// =======================

class DataLine2D{ public:
    // data
    int      n = 0;
    bool   bSharedX =false;
    //Arr xs;
    double * xs    =0;
    double * ys    =0;
    Func1d   yfunc =0;
    //update_data    =false;
    // draw properties
    Rect2d   bounds;
    int      glObj      =0;
    char     lineStyle  ='-';
    //char     pointStyle ='+';
    char     pointStyle =' ';
    float    pointSize  =0.1;
    bool     bView = true;
    uint32_t clr        =0xFFFF00FF;
    std::string label;


    void update();
    int  render();
    void draw();
    void view();

    inline void set     ( int n_, double * xs_, double * ys_){ n=n_; xs=xs_; ys=ys_; };
    //inline void allocate( int n_ ){ n=n_; xs=new double[n]; ys=new double[n]; bSharedX=false; };
    inline void allocate( int n_ ){ n=n_; _allocIfNull(xs,n); _allocIfNull(ys,n); bSharedX=false; };

    inline void linspan(double xmin, double xmax){ VecN::linspan(n,xmin,xmax,xs); };
    inline void arange (double xmin, double dx  ){ VecN::arange (n,xmin,dx,xs);   };

    inline DataLine2D()=default;
    inline DataLine2D(int n_){ allocate(n_); }
    inline DataLine2D(int n_,double*xs_){ n=n_; bSharedX=true; xs=xs_; ys=new double[n]; }
    //inline DataLine2D(int n_,double xmin,double xmax){ allocate(n_); linspan(xmin,xmax);  }
    inline DataLine2D(int n_,double xmin,double dx, uint32_t clr_=0xFFFF00FF, std::string label_="", double* ys_=0 ){ ys=ys_; allocate(n_); arange(xmin,dx);                 clr=clr_; label=label_; }
    inline DataLine2D(int n_,double*xs_,            uint32_t clr_=0xFFFF00FF, std::string label_="", double* ys_=0 ){ ys=ys_; n=n_; bSharedX=true; xs=xs_; ys=new double[n]; clr=clr_; label=label_; }

    ~DataLine2D();
};

template<typename Func>
void evalLine( DataLine2D& l, Func func ){
    for(int i=0; i<l.n; i++){ l.ys[i]=func(l.xs[i]); }
};

// =======================
// ====  Plot2D
// =======================

class Plot2D{ public:
    // data
    std::vector<DataLine2D*> lines;
    // properties
    Rect2d   bounds;
    int      nXTicks=0,nYTicks=0;
    double * xTicks=NULL;
    double * yTicks=NULL;
    int      glObj=0;


    Vec2d shift   =Vec2dZero;
    Vec2d scaling =Vec2dOnes;

    Vec2d  axPos;
    Rect2d axBounds;
    float tickSz = 0.1;
    Vec2f legend_pos={0.,0.};

    std::string xlabel="";
    std::string ylabel="";

    inline DataLine2D* add(DataLine2D* dline){
        lines.push_back(dline);
        return dline;
    };

    bool     logX  = false;
    bool     logY  = false;
    bool     bGrid       = true;
    bool     bAxes       = true;
    bool     bTicks      = true;
    bool     tickCaption = false;
    uint32_t clrBg      = 0x00f0f0f0;
    uint32_t clrGrid    = 0xFF858585;
    uint32_t clrTicksX  = 0xFF000000;
    uint32_t clrTicksY  = 0xFF000000;
    char   * tickFormat = "%2.2f\0";
    int      fontTex=0;



    void update();
    void drawTexts();
    void drawAxes();
    int  render();
    void view  (bool bAxes=true);
    void init  ();
    //void xsharingLines(int nl, int np);
    //void xsharingLines(int nl, int np, double xmin, double dx);
    //void init( double dx, double dy );
    void xsharingLines(int nl, int np, double xmin=0.0, double dx=0.1, uint32_t* colors=0, int ncol=-1);
    void autoAxes(double dx, double dy);
    void clear(bool bDeep=true);
    void erase(int i);

    void savetxt(const char* fname);

    void drawHline ( double y );
    void drawVline ( double x );
    //void drawCursor( Vec2d p, double sz );

    //void addDataLineSequence( int nl, int np, uint32_t colors, int ncol=-1 ){}

    Plot2D()=default;
    ~Plot2D(){};
};

// =======================
// ====  QuePlot2D
// =======================

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

