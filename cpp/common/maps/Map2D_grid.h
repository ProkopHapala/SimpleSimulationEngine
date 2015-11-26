class Map2D_grid{
	public:
	int pow;
	int pow_nt;
	int pow_tn;

	int n;
	int nt;
	int tn;

	int nt2;
	int tn2;

	MapTile2D ** tiles;

	Map2D_grid( int pow_nt_, int pow_tn_ ){
		int pow_nt = pow_tn_;
		int pow_tn = pow_nt_;
		int pow    = pow_nt + pow_tn;
		n  = 1 << pow;
		nt = 1 << pow_nt;
		tn = 1 << pow_tn;
		nt2 = nt*nt;
		tn2 = tn*tn;
		tiles = new MapTile2D*[ nt2 ];
		for(int it=0; it<nt2; it++){
			tiles[it] = NULL;
		}
	}

	inline int getTileIndex( int ix, int iy ){ 
		int itx  = ix >> pow_tn;   
		int ity  = iy >> pow_tn;   
		return ( ity << pow_nt ) + itx;		
    }

	inline void getIndexes( int ix, int iy, int& iTile, int& iCell ){
		int itx  = ix >> pow_tn;   
		int ity  = iy >> pow_tn;   
		iTile = ( ity << pow_nt ) + itx;	
		int tix = ix - ( itx << pow_tn ); 
		int tiy = iy - ( ity << pow_tn );
		iCell = tiy << pow_tn + tix;
	}; 

	inline     MapCell2D * getCell( int ix, int iy ){ 
		int iTile,iCell;
		getIndexes( ix, iy, iTile, iCell );
		return &tiles[ iTile ]->cells[ iCell ];  
	};

};
