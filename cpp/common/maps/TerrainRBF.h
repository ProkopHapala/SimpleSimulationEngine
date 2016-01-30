
#ifndef  TerrainRBF_h
#define  TerrainRBF_h

//#include "fastmath.h"
//#include "Vec2.h"

#include "RBF2D.h"
#include "HashMap2D.h"

class TerrainRBF : public HashMap2D<RBF2D> {
    public:
    bool exclusive_rbfs = false; // this is true if rbfs shoudl be deleted with this terrain
    int nrbfs;
    RBF2D * rbfs = NULL;

    int nbuff;
    RBF2D ** buff = NULL;

    // ==== functions
    double getVal        ( double x, double y );
	int    renderRect    ( double xmin, double ymin, double xmax, double ymax, int nx );
    int    insertRBF     ( RBF2D* rbf );
    int    removeRBF     ( RBF2D* rbf );
    int    insertRBFs    ( int n_, RBF2D * rbfs_ );
    void   generateRandom( double xmin, double ymin, double xmax, double ymax, double rmin, double rmax, double vmin, double vmax, int n_ );

    // ==== inline functions

    void init( double step_, UINT power_, int nbuff_ ){
        HashMap2D<RBF2D>::init( step_, power_ );
        nbuff = nbuff_;
        buff  = new RBF2D*[nbuff];
    }

    ~TerrainRBF(){
        if ( ( rbfs != NULL ) && exclusive_rbfs ) delete rbfs;
        if (   buff != NULL )                     delete buff;
    }

};

#endif

