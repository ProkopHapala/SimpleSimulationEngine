#ifndef  TileTree2D_h
#define  TileTree2D_h

/*
 TO DO
   - is there any elegant way how to make it general ( i.e. for any order of branching; not just 2 levels ) without loss of performance ?
   - how to implement removel of tile ? Would depend on wheather OBJECT is pointer or not;
*/

//template < class OBJECT, unsigned int POWER >
template < class OBJECT, int POWER >
class LeafTile2D{
	public:
    constexpr unsigned static int n    = 1<<POWER;
	constexpr unsigned static int n2   = n*n;

    //constexpr static int n    = 1<<POWER;
	//constexpr static int n2   = n*n;

	int nfilled; // this is used mostly for removing

	OBJECT cells[n2];

	// ==== inline functions

	inline int index2D       ( int ix, int iy ) const { return ( iy << POWER ) + ix;      };
	inline OBJECT* getPointer( int ix, int iy )       { return &cells[ index2D(ix,iy) ];  };

	inline void fill( OBJECT obj ){ for( int i=0; i<n2; i++ ){ cells[i]=obj; } }

	LeafTile2D(){ nfilled = 0; }

};

template < class OBJECT, unsigned int POWER, unsigned int NX, unsigned int NY >
//template < class OBJECT, int POWER, int NX, int NY >
class TileTree2D{
	public:
/*
    constexpr unsigned static int power  = POWER;
    constexpr unsigned static int nx     = NX;
    constexpr unsigned static int ny     = NY;

    constexpr unsigned static int nsub      = 1<<POWER;
	constexpr unsigned static int sub_mask  = nsub - 1;
	constexpr unsigned static int nsub2     = nsub*nsub;
	constexpr unsigned static int nxy       = NX*NY;

	constexpr unsigned static int ntotx  = nsub  * NX;
	constexpr unsigned static int ntoty  = nsub  * NY;
	constexpr unsigned static int ntotxy = ntotx * ntoty;
*/


    constexpr static int power  = POWER;
    constexpr static int nx     = NX;
    constexpr static int ny     = NY;

    constexpr static int nsub      = 1<<POWER;
	constexpr static int sub_mask  = nsub - 1;
	constexpr static int nsub2     = nsub*nsub;
	constexpr static int nxy       = NX*NY;

	constexpr static int ntotx  = nsub  * NX;
	constexpr static int ntoty  = nsub  * NY;
	constexpr static int ntotxy = ntotx * ntoty;

	LeafTile2D<OBJECT,POWER>* tiles[ nxy ];

	// ==== inline functions

	inline int getSubIndex( int i          ) const { return   i  & sub_mask;   };
	inline int getSupIndex( int i          ) const { return   i  >> POWER;     };
	inline int isup2D     ( int ix, int iy ) const { return ( iy *  NX ) + ix; };

	inline OBJECT* getPointer( int ix, int iy ) {
		int i           = isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		LeafTile2D<OBJECT,POWER>* tile = tiles[ i ];
		if ( tile == NULL ) return NULL;
		return tile->getPointer( getSubIndex( ix ), getSubIndex( iy ) );
	}

	inline OBJECT* getValidPointer( int ix, int iy, OBJECT null_val ){
		int i           = isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		LeafTile2D<OBJECT,POWER>* tile = tiles[ i ];
		if ( tile == NULL ){
			tile = new LeafTile2D<OBJECT,POWER>();
			tiles[ i ] = tile;
			tile->fill( null_val );
		}
		return tile->getPointer( getSubIndex( ix ), getSubIndex( iy ) );
	}

	TileTree2D( ){
		for( int i=0; i<nxy; i++ ){ tiles[i] = NULL; };
	}

	~TileTree2D(){
		for( int i=0; i<nxy; i++ ){ if( tiles[i] != NULL ) delete tiles[i]; };
	}

};


template < class OBJECT, unsigned int POWER, unsigned int NX, unsigned int NY >
class TileTree2D_d : public TileTree2D<OBJECT,POWER,NX,NY> {

    using PARENT = TileTree2D<OBJECT,POWER,NX,NY>;

    public:
	double xmin,ymin,xmax,ymax,xspan,yspan;
	double xStep,yStep;
	double invXStep,invYStep;

	inline int    getIx ( double  x ) const{ return (int)( invXStep * ( x - xmin ) ); };
	inline int    getIy ( double  y ) const{ return (int)( invYStep * ( y - ymin ) ); };
	inline double getX  ( int    ix ) const{ return ( xStep * ix ) + xmin;                };
	inline double getY  ( int    iy ) const{ return ( yStep * iy ) + ymin;                };

	inline OBJECT* getPointer_d     ( double x, double y                  ){ return PARENT::getPointer( getIx(x), getIy(y)           );       };
	inline OBJECT* getValidPointer_d( double x, double y, OBJECT null_val ){ return PARENT::getValidPointer( getIx(x), getIy(y), null_val );  };

	void setPos( double xmin_, double ymin_ ){
        xmin = xmin_; ymin = ymin_;
        xmax  = xspan + xmin;
		ymax  = yspan + ymin;
	}

    void setSpacing( double yStep_, double xStep_ ){
		xStep = xStep_; invXStep = 1/xStep;
		yStep = yStep_; invYStep = 1/yStep;
		xspan = xStep * PARENT::ntotx;
		yspan = yStep * PARENT::ntoty;
        setPos( -xspan * 0.5d, -yspan * 0.5d );
	}

};


#endif
