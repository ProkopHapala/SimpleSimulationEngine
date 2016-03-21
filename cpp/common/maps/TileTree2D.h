#ifndef  TiledView_h
#define  TiledView_h

/*
 TO DO  
   - is there any elegant way how to make it general ( i.e. for any order of branching; not just 2 levels ) without loss of performance ?
   - how to implement removel of tile ? Would depend on wheather OBJECT is pointer or not;

*/

template < class OBJECT, unsigned int POWER >
class LeafTile2D{
	public:
	//int power;
    constexpr unsigned static int n    = 1<<POWER; 
	//constexpr unsigned static int mask = n - 1;
	constexpr unsigned static int n2   = n*n;

	int nfilled; // this is used mostly for removing

	OBJECT cells[n2];

	// ==== inline functions

	inline int index2D       ( int ix, int iy ){ return ( iy << POWER ) + ix;      };
	inline OBJECT* getPointer( int ix, int iy ){ return &cells[ index2D(ix,iy) ];  };

/*
	inline void getBare      ( int ix, int iy,       OBJECT& obj ){ obj = cells[ index2d(ix,iy) ];  }; 
	inline void setBare      ( int ix, int iy, const OBJECT& obj ){ cells[ index2d(ix,iy) ] = obj;  };

	inline bool get       ( int ix, int iy,       OBJECT& obj ){ 
		obj = cells[ index2d(ix,iy) ];  
		return obj.isEmpty();
	}; 

	inline void set       ( int ix, int iy, const OBJECT& obj ){ 
		int  i = index2d(ix,iy);
		char emptyness = ( cells[ i ].isEmpty()  || obj.isEmpty()<<1 );
		if       ( emptyness == 1 ){
			nfilled++;
		}else if ( emptyness == 2 ){
			nfilled--;
		}
		cells[ i ] = obj;              
	};
*/

	inline void fill( OBJECT obj ){ for( int i=0; i<n2; i++ ){ cells[i]=obj; } }

	LeafTile2D(){ nfilled = 0; }

};

/*
template <class TILE, unsigned int POWER >
class BranchTile2D{
	public:
	//int power;
    constexpr unsigned static int n    = 1<<POWER; 
	constexpr unsigned static int mask = n - 1;
	constexpr unsigned static int n2   = n*n;

	int nfilled; // this is used mostly for removing

	TILE cells[n2];

	// ==== inline functions

	inline int index2d    ( int ix, int iy                    ){ return ( iy << POWER ) + ix;    };
	inline void getBare   ( int ix, int iy,       OBJECT& obj ){ obj = cells[ index2d(ix,iy) ];  }; 
	inline void setBare   ( int ix, int iy, const OBJECT& obj ){ cells[ index2d(ix,iy) ] = obj;  };

	inline bool get       ( int ix, int iy,       OBJECT& obj ){ 
		obj = cells[ index2d(ix,iy) ];  
		return obj.isEmpty();
	}; 

	inline void set       ( int ix, int iy, const OBJECT& obj ){ 
		int  i = index2d(ix,iy);
		char emptyness = ( cells[ i ].isEmpty()  || obj.isEmpty()<<1 );
		if       ( emptyness == 1 ){
			nfilled++;
		}else if ( emptyness == 2 ){
			nfilled--;
		}
		cells[ i ] = obj;              
	};

	LeafTile2D(){ nfilled = 0; }

}
*/

//template < template TILE<OBJECT>, class OBJECT, unsigned int POWER, unsigned int NX, unsigned int NY >
template < class OBJECT, unsigned int POWER, unsigned int NX, unsigned int NY >
class TileTree2D{
	public:
    constexpr unsigned static int nsub      = 1<<POWER; 
	constexpr unsigned static int sub_mask  = nsub - 1;
	constexpr unsigned static int nsub2     = nsub*nsub;
	constexpr unsigned static int nxy       = NX*NY;

	constexpr unsigned static int ntotx  = nsub  * NX;
	constexpr unsigned static int ntoty  = nsub  * NY;
	constexpr unsigned static int ntotxy = ntotx * ntoty;

	LeafTile2D<OBJECT,POWER>* tiles[ nxy ];

	// ==== inline functions

	inline int getSubIndex( int i          ){ return   i  & sub_mask;  };
	inline int getSupIndex( int i          ){ return   i  >> POWER;   };
	inline int isup2D     ( int ix, int iy ){ return ( iy *  NX ) + ix; };

	inline OBJECT* getPointer( int ix, int iy ){
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

/*
	inline bool getBare( int ix, int iy,       OBJECT& obj ){
		int i           = isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] != NULL ){	return tiles[ i ]->getBare( getSubIndex( ix ), getSubIndex( iy ), obj );	}
	}

	inline bool setBare( int ix, int iy, const OBJECT& obj ){
		int i= isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] = NULL ){ tiles[ i ] = new LeafTile2D<OBJECT,POWER>(); }
		tiles[ i ]->setBare( getSubIndex( ix ), getSubIndex( iy ) );
	}

	inline bool get( int ix, int iy,       OBJECT& obj ){
		int i           = isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] != NULL ){
			return tiles[ i ]->get( getSubIndex( ix ), getSubIndex( iy ), obj );
		}
		return false;
	}

	inline bool set( int ix, int iy, const OBJECT& obj ){
		int i= isup2D( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] = NULL ){ tiles[ i ] = new LeafTile2D<OBJECT,POWER>(); }
		tiles[ i ]->set( getSubIndex( ix ), getSubIndex( iy ) );
		if( tiles[ i ]->nfilled == 0 ){ delete tiles[ i ]; }
	}
*/

	TileTree2D( ){
		//tiles = LeafTile2D<OBJECT,POWER>( );
		for( int i=0; i<nxy; i++ ){ tiles[i] = NULL; };
	}

	~TileTree2D(){
		for( int i=0; i<nxy; i++ ){ if( tiles[i] != NULL ) delete tiles[i]; };
		//delete tiles;
	}

};


#endif
