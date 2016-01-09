
#ifndef  GridMap2D_h
#define  GridMap2D_h

#include <forward_list>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "MapSite2D.h"

// ========== GridMap

template <class TYPE1> 
class GridMap2D {
	public:
	double step, invStep;
	int nx, ny, nxy;
	bool isStatic = false;
	MapSite2D * store;

	// ==== function declarations

	void makeStatic( );
	int  insertTriangle( Triangle2D* t );
	int  insertLine( Line2D* l );

	// ==== inline functions
	
	inline int getIx( double x ){ return (int)( invStep * x ); };
	inline int getIy( double y ){ return (int)( invStep * y ); };

	inline int getIndex( int ix, int iy     ){ return nx*iy + ix; }
	inline int getIndex( double x, double y ){ return getIndex( getIx(x), getIy(y) ); };

	inline void insert( TYPE1* p, int i ){
		store  [ i ] -> push_front( p );
		ns     [ i ] ++;
	};
	inline void insert( TYPE1* p, int ix, int iy ){ insert( p, getIndex( ix,iy ) ); };
	
};


template <class TYPE1> 
void GridMap2D::init( int nx_, int ny_, double step_ ){
	step    = step_; 
	invStep = 1/step;
	nx = nx_; ny=ny_;
	nxy = nx*ny;
	ns     = new int[ nxy ];
	store  = new std::forward_list<TYPE1*>*[ nxy ];
	for (int i=0; i<nxy; i++){ store[i] = new std::forward_list<TYPE1*>(); ns[i] = 0; }
};

template <class TYPE1> 
void GridMap2D::makeStatic( ){
	store_static = new TYPE1**[nxy];
	for (int i=0; i<nxy; i++){ 
		int ni = ns[i];
		if ( ni>0 ){ 
			//printf( " ni %i \n", ni );
			TYPE1** arr = new TYPE1*[ ni ];
			store_static[ i ] = arr;
			auto list = store [ i ];
			int j = 0;
			for ( auto it = list->cbegin(); it != list->cend(); ++it){
				arr[ j ] = *it;
				j++;
			}  
		}
	}
	isStatic = true;
	// dealocation
	for (int i=0; i<nxy; i++){ delete store[i]; }
	delete [] store;
}

#endif 

