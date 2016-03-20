#ifndef  TiledView_h
#define  TiledView_h

/*
 TO DO  
   - is there any elegant way how to make it general ( i.e. for any order of branching; not just 2 levels ) without loss of performance ?
   - how to implement removel of tile ? Would depend on wheather OBJECT is pointer or not;

*/

template <class OBJECT, unsigned int POWER >
class LeafTile2D{
	public:
	//int power;
    constexpr unsigned int n    = 1<<POWER; 
	constexpr unsigned int mask = n - 1;
	constexpr unsigned int n2   = n*n;

	int nfilled; // this is used mostly for removing

	OBJECT cells[n2];

	inline int index2d( int ix, int iy ){ return ( iy << POWER ) + ix; };

	inline bool get( int ix, int iy,       OBJECT& obj ){ obj = cells[ index2d(ix,iy) ]; return true; } // we should inform if tile is empty
	inline void set( int ix, int iy, const OBJECT& obj ){ cells[ index2d(ix,iy) ] = obj;              }

}

template <TILE,OBJECT,unsigned int POWER>
class BranchTile2D{
	public:
	int sub_pow;
	int sub_mask;
	int nsub,nsub2;

	int ntotx,ntoty,ntotxy;
	int nx,ny,nxy;

	TILE<OBJECT> * tiles;

	inline int getSubIndex( int i          ){ return   i   & sub_mask;  };
	inline int getSupIndex( int i          ){ return   i  >> sub_pow;   };
	inline int isup2D     ( int ix, int iy ){ return ( iy *  nx ) + ix; };

	inline bool get( int ix, int iy,       OBJECT& obj ){
		int i           = index2d( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] != NULL ){
			return tiles[ i ].get( getSubIndex( ix ), getSubIndex( iy ), obj );
		}
		return false;
	}

	inline bool set( int ix, int iy, const OBJECT& obj ){
		int i= index2d( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] = NULL ){ tiles[ i ] = new TILE<OBJECT>( sub_pow ); }
		tiles[ i ].set( getSubIndex( ix ), getSubIndex( iy ) );
		if( tiles[ i ].nfilled == 0 ){ delete tiles[ i ]; }
	}

	BranchTile2D{

	}

	~BranchTile2D(){
		for( int i=0; i<nxy; i++ ){ delete tiles };
		delete tiles;
	}

}


#endif
