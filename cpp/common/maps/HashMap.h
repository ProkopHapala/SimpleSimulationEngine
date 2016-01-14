
#ifndef  HashMap_h
#define  HashMap_h

typedef unsigned long    ulong;
typedef unsigned int     uint;


// TODO: maybe optimal would be put several objects per bucket

template <class TYPE1 > 
HashMapField{
	public:
	TYPE1*  object;  // pointer is also unique identifier
	ulong   bucket;  // unique index of bucket 
	uint    n;       // number of stored objects with this hash 
	//uint    hash;    // hash    // not necessary to store, can be quickly computed from ibox

	HashMapField(){
		object = NULL;
		n      = 0; 
	}

	inline void set( TYPE1* object_, ulong bucket_ ){
		object  = object_;
		bucket  = bucket_;
		//hash    = hash_;
	}

}

template <class TYPE > 
class HashMap{
	public:
	uint power;
	uint mask;
	uint capacity;
	uint filled=0;

	int DEBUG_counter;

	HashMapField<TYPE>*  fields;

	inline ulong  hashFunc( unsigned long index ){ 
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

	void init( uint power_ ){
		power      = power_;
		capacity   = 1<<power;
		mask       = capacity-1;
		//printf( "  power %i capacity %i mask %i \n ", power, capacity, mask );
		fields     = new HashMapField<TYPE>[ capacity ];
		// initialization done by constructor (?)
	}

	inline void insert( TYPE* object, ulong bucket, uint hash, uint i ){
		fields[ hash ].n--;
		fields[ i    ].set( object, bucket ); 
		filled++;
	} 
	inline uint insert( TYPE* object, ulong bucket ){
		if ( (filled<<1) > capacity ) resize( power+1 );
		uint hash = hashFunc( bucket );   // must be after resize
		uint i = hash;
		while( fields[i].object != NULL ) i=(i+1)&mask; 
		insert( object, bucket, hash, i );
		return i;
	};

	void remove( uint i, uint hash ){
		filled--;
		fields[ hash ].n--;
		fields[ i    ].object = NULL; 
	};

	void resize( uint power_ ){
		uint old_capacity = capacity;
		HashMapField<TYPE>*  old_fields = fields;
		init( power_ );
		HashMapField<TYPE>*  p_field    = old_fields;
		for (int i=0; i<old_capacity; i++){ 
			insert( p_field->object, p_field->bucket );
			p_field++;
		}
		delete old_fields;
	}

	int find( TYPE* object, ulong bucket ){
		uint hash = hashFunc( bucket );
		uint n    = fields[ hash ].n;
		uint i    = hash;
		while( n > 0 ){
			if( object == fields[ i ].object             ){ return i; }
			if( hash   == hashFunc( fields[ i ].bucket ) ){ n--; } // FIXME : here we call hashFunc in loop, would storing hash improve performance ?  
			i=(i+1)&mask;
		}
		return -1;
	}

	int getBucket( ulong bucket, uint * outi ){
		uint hash = hashFunc( bucket );
		uint n    = fields[ hash ].n;
		uint i    = hash;
		uint j    = 0;
		while( n > 0 ){
			ulong bucketi = fields[ i ].bucket;
			if( hash == hashFunc( bucketi ) ){   // FIXME : here we call hashFunc in loop, would storing hash improve performance ? 
				if( bucketi == bucket ){
					outi[ j ] = i;
					j++;
				}
				n--;
			}   
			i=(i+1)&mask;
		}
		return -1;
	}

};


#endif

