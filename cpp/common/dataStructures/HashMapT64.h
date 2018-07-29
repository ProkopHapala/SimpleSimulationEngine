#ifndef HashMapT64_h
#define HashMapT64_h

/*
NOTES:

Size:
 - 64 bit boxId is 21 bit/axis = 2,097,152;   with 1km per box this is just 2 milion km per axis, which does not even cover Jupiter system
 - 32 bit per axis will allow 4 billion ticks ( 4 billion km per axis with 1km resolution ) that is approx distance to pluto
 - 64 bit per axis will allow 18.446744 billion kilometers (1.949 ly ) with 1 mm resolution  (light year 9.4607e+15 m)
 - 52 bit mantisa of double means 4.5e+15 ticks ( 0.5 ly with 1 m resolution ), distance to Pluto (4.43682e+9 km) with 1mm resolution
 - anyway trepasing >1000 boxes per raycast is super inefficient

*/

template<typename T>
class HashMapT64{ public:
	uint8_t  power    =0;
	uint32_t capacity =0;
	uint32_t mask     =0;
	uint32_t filled   =0;
	float    maxFillRatio = 0.6;
	float    begFillRatio = 0.3;

    T        EMPTY_O; // TO BE SET
	uint32_t EMPTY_H = 0xFFFFFFFF;

	int DEBUG_counter=0;

	T*         store;  // pointer or index of stored object
	uint32_t*  hashs;  // hash   // not necessary to store, can be quickly computed from ibox
	uint16_t*  hits;   // number of objects with this hash
	//uint8_t*  hits;  // number of objects with this hash

	inline uint32_t hash( uint64_t ibox ){
		// see http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
		// https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key
		//return  ( 2166136261UL ^ i * 16777619 ) & mask;
		//return  ( 2166136261UL ^ (i * 16777619) ) & mask;
		//i = 1664525*(i+15445) ^ 1013904223; return (1664525*i ^ 1013904223) & mask;
		//return ((i >> power)^i) & mask;
		//return (i*2654435761 >> 16)&mask;   // Knuth's multiplicative method
		return ( ibox * 11400714819323198549L )&mask;  // Knuth's multiplicative method  64bit
		//i = ((i >> 16) ^ i) * 0x45d9f3b; i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x3335b369; i = ((i >> 16) ^ i) * 0x3335b369; return  ((i >> 16) ^ i)&mask;
	}

	inline void set( int i, const T& p, uint32_t h ){
		store[ i ] =  p;
		hashs[ i ] =  h;
	}

	void init( int power_, T EMPTY_O_ ){
        EMPTY_O    = EMPTY_O;
		power      = power_;
		capacity   = 1<<power;
		mask       = capacity-1;
		//printf( "  power %i capacity %i mask %i \n ", power, capacity, mask );
		store      = new T   [capacity];
		hits       = new int [capacity];
		hashs      = new int [capacity];
		for (int i=0; i<capacity; i++){
			hits  [i] =  0;
			set( i, EMPTY_O, EMPTY_H );
		}
	}

	inline uint32_t getAllInBox( uint64_t ibox, uint32_t h, int* outi ){
		DEBUG_counter =0;
		uint32_t n = hits[h];
		uint32_t i = h;
		uint32_t j = 0;
		while( n>0 ){
			if( hashs[i]==h ){
				if( store[i].getCellIndex() == ibox ){
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
	inline uint32_t getAllInBox( uint64_t ibox, int* outi ){ uint32_t h = hash( ibox ); return getAllInBox( ibox, h, outi ); }

	inline uint32_t insert( const T& o, uint32_t h, uint32_t i ){
		hits   [h] ++;
		set( i, o, h );
		filled++;
	}
	inline uint32_t insert( const T& o, uint32_t h ){
		checkResize();
		uint32_t i = h;
		while( !store[i].isEmptyCell() ) i=(i+1)&mask;
		insert( o, h, i );
		return i;
	};
	inline uint32_t insert( const T& o ){ uint32_t h = hash( o.getCellIndex() ); return insert( o, h );   }

	void remove( int i, uint32_t h ){
		filled--;
		hits  [ hash[i] ] -- ;
		set( i, EMPTY_O, EMPTY_H );
	};

	void resize( int power_ ){
		uint32_t   old_capacity = capacity;
		uint32_t*  old_hashs    = hashs;
        T*         old_store    = store;
		delete [] hits;
		init( power_ );
		for (int i=0; i<old_capacity; i++){
			insert( old_store[i], old_hashs[i] );
		}
		delete [] old_store;
		delete [] old_hashs;
	}

	void checkResize(){
	   if( filled > capacity * maxFillRatio ){
	        power++;
            while( (1<<power)*begFillRatio < filled ) power++;
	        resize( power );
	   };
	}

	void reserve( int n ){
        int nTarget = (int)(n/begFillRatio)+1;
        if( nTarget > capacity ){
            while( (1<<power) < nTarget )power++;
            resize( power );
        }
	}

};

#endif

