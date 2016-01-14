

template <class TYPE1 > class HashMap{
	public:
	int power;
	int mask;
	int capacity;
	int filled=0;

	int DEBUG_counter;

	TYPE1*  store;
	int*    hits;  // number of objects with this hash 
	int*    iboxs; // unique box index 
	int*    hashs; // hash   // not necessary to store, can be quickly computed from ibox

	inline int  hash( unsigned int i ){ 
		// see http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key 
		//return  ( 2166136261UL ^ i * 16777619 ) & mask; 
		//return  ( 2166136261UL ^ (i * 16777619) ) & mask; 
		//i = 1664525*(i+15445) ^ 1013904223; return (1664525*i ^ 1013904223) & mask;
		//return ((i >> power)^i) & mask;
		return (i*2654435761 >> 16)&mask;   // Knuth's multiplicative method
		//i = ((i >> 16) ^ i) * 0x45d9f3b; i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x3335b369; i = ((i >> 16) ^ i) * 0x3335b369; return  ((i >> 16) ^ i)&mask;
	}

	inline void set( int i, TYPE1 p, int ibox, int h ){
		store  [ i ] =  p; 
		iboxs  [ i ] =  ibox;
		hashs  [ i ] =  h;
	}

	void init(	int power_ ){
		power      = power_;
		capacity   = 1<<power;
		mask       = capacity-1;
		//printf( "  power %i capacity %i mask %i \n ", power, capacity, mask );
		store      = new TYPE1 [capacity]; 
		hits       = new int   [capacity]; 
		hashs      = new int   [capacity];
		iboxs      = new int   [capacity];
		for (int i=0; i<capacity; i++){ 
			hits  [i] =  0    ; 
			set( i, NULL, 0xFFFFFFFF, 0xFFFFFFFF ); 
		}
	}

	inline int getAllInBox( int ibox, int h, int* outi ){
		DEBUG_counter =0;
		int n = hits[h];
		int i = h;
		int j =0;
		while( n>0 ){
			if( hashs[i]==h ){ 
				if( iboxs[i]==ibox ){
					outi[j] = i; 
					j++;
				}
				n--; 
			};
			i=(i+1)&mask;
			DEBUG_counter++;
		}
		return j;
	}
	inline int getAllInBox( int ibox, int* outi ){ int h = hash( ibox ); return getAllInBox( ibox, h, outi ); }

	inline int getAllInBox_noHash( int ibox, int* outi ){
		DEBUG_counter =0;
		int h = hash( ibox );
		int n = hits[h];
		int i = h;
		int j =0;
		while( n>0 ){
			int ib = iboxs[i];
			int ih = hash( ib );
			if( ih==h ){ 
				if( ib==ibox ){
					outi[j] = i; 
					j++;
				}
				n--; 
			};
			i=(i+1)&mask;
			DEBUG_counter++;
		}
		return j;
	}

	inline int insert( TYPE1 p, int ibox, int h, int i ){
		hits   [h] ++;
		set( i, p, ibox, h ); 
		filled++;
	} 
	inline int insert( TYPE1 p, int ibox, int h ){
		if ( (filled<<1) > capacity ) resize( power+1 );
		int i = h;
		while( store[i] != NULL  ) i=(i+1)&mask; 
		insert( p, ibox, h, i );
		return i;
	};
	inline int insert( TYPE1 p, int ibox ){ int h = hash( ibox ); return  insert( p, ibox, h );   }

	void remove( int i, int h ){
		filled--;
		hits  [ hash[i] ] -- ;
		set( i, NULL, 0xFFFFFFFF, 0xFFFFFFFF ); 
	};

	void resize( int power_ ){
		int old_capacity = capacity;
		TYPE1*  old_store = store; 
		int*    old_hashs = hashs;
		int*    old_iboxs = iboxs;
		delete hits;
		init( power_ );
		for (int i=0; i<old_capacity; i++){ 
			insert( old_store[i], old_iboxs[i], old_hashs[i] );
		}
		delete old_store;
		delete old_hashs;
		delete old_iboxs;
	}

	void writeToFile( FILE * pFile ){
		char str[ 1000 ];
		for ( int i=0 ; i<capacity ; i++)   {
			if( store[i] != NULL ){ 
				store[i]->toString( str ); 				
				fprintf ( pFile, "%s \n", str ); 
			}
		}
	}

};

