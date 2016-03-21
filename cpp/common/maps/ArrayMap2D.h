#ifndef  ArrayMap2D_h
#define  ArrayMap2D_h

template < class OBJECT, int NX, int NY >
class ArrayMap2D{
	public:

    constexpr static int nx     = NX;
    constexpr static int ny     = NY;
	constexpr static int nxy    = NX*NY;

	OBJECT cells[ nxy ];

    public:
	double xmin,ymin,xmax,ymax,xspan,yspan;
	double xStep,yStep;
	double invXStep,invYStep;

	// ==== inline functions INDEPENDENT of real-space sampling

	inline int isup2D ( int ix, int iy ) const { return ( iy *  NX ) + ix; };

	// ==== inline functions DEPENDENT of real-space sampling

	inline int    getIx ( double  x ) const{ return (int)( invXStep * ( x - xmin ) ); };
	inline int    getIy ( double  y ) const{ return (int)( invYStep * ( y - ymin ) ); };
	inline double getX  ( int    ix ) const{ return ( xStep * ix ) + xmin;            };
	inline double getY  ( int    iy ) const{ return ( yStep * iy ) + ymin;            };

	inline OBJECT* getPointer ( double x, double y  ){ return &cells[ isup2D( getIx(x), getIy(y) ) ]; };

	void setPos( double xmin_, double ymin_ ){
        xmin = xmin_; ymin = ymin_;
        xmax  = xspan + xmin;
		ymax  = yspan + ymin;
	}

    void setSpacing( double yStep_, double xStep_ ){
		xStep = xStep_; invXStep = 1/xStep;
		yStep = yStep_; invYStep = 1/yStep;
		xspan = xStep * nx;
		yspan = yStep * ny;
        setPos( -xspan * 0.5d, -yspan * 0.5d );
	}

};

#endif
