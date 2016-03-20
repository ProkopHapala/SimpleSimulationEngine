#ifndef  TiledView_h
#define  TiledView_h

/*
 TO DO  
   - is there any elegant way how to make it general ( i.e. for any order of branching; not just 2 levels ) without loss of performance ?
   - how to implement removel of tile ? Would depend on wheather OBJECT is pointer or not;

*/

template <class OBJECT >
class LeafTile2D{
	public:
	int power;
	int mask;
	int n,n2;

	OBJECT * cells;

	inline int index2d    ( int ix, int iy ){ return ( iy << power ) + ix; };

	inline OBJECT get( int ix, int iy         ){ return cells[ index2d(ix,iy) ];   }
	inline void   set( int ix, int iy , OBJECT){ cells[ index2d(ix,iy) ] = OBJECT; }

	OBJECT * tiles;

}

template <OBJECT>
class BranchTile2D{
	public:
	int sub_pow;
	int sub_mask;
	int nsub,nsub2;

	int ntotx,ntoty,ntotxy;
	int nx,ny,nxy;

	LeafTile2D * tiles;

	inline int getSubIndex( int i          ){ return   i   & sub_mask;  };
	inline int getSupIndex( int i          ){ return   i  >> sub_pow;   };
	inline int isup2D     ( int ix, int iy ){ return ( iy *  nx ) + ix; };

	inline TILE get( int ix, int iy ){
		int i           = index2d( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] != NULL ){
			return tiles[ i ].get( getSubIndex( ix ), getSubIndex( iy ) );
		}
	}

	inline bool set( int ix, int iy, TILE ){
		int i= index2d( getSupIndex( ix ), getSupIndex( iy ) );
		if( tiles[ i ] = NULL ){ ;                           }
		tiles[ i ].set( getSubIndex( ix ), getSubIndex( iy ) );
	}

}


#endif
