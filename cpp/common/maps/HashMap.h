
#ifndef  HashMap_h
#define  HashMap_h

#include <stdio.h>


//typedef unsigned long    ULONG;
typedef unsigned int     ULONG;
typedef unsigned int     UINT;

inline ULONG  hashFunc( ULONG i ) {
	// see http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
	//return  ( ( 2166136261UL ^ i ) * 16777619 );  // VERY BAD
	//return  ( 2166136261UL ^ (i * 16777619) );    // VERY BAD
	//return ((i >> 16)^i);                           // NOT WORKING
	//return ( i * 2654435761 >> 16 );   // Knuth's multiplicative method // GOOD
	return ( i * 2654435761 >> 8 );   // GOOD, even better than Knuth for map size 2^16 ( 65536 fields )
	//return ( i * 2654435761 >> 16 ) ^ i; // REASONABLY GOOD but worse than pure Knuth
	//i = ((i >> 16) ^ i) * 0x45d9f3b; i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i); // REASONABLY GOOD but worse than pure Knuth
	//i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i); // LESS GOOD
	//i = ((i >> 16) ^ i) * 0x3335b369; i = ((i >> 16) ^ i) * 0x3335b369; return  ((i >> 16) ^ i); // LESS GOOD
}

// TODO: maybe optimal would be put several objects per bucket

template <class TYPE1 >
class HashMapField{
	public:
	TYPE1*  object;  // pointer is also unique identifier
	ULONG   bucket;  // unique index of bucket
	UINT    n;       // number of stored objects with this hash
	//UINT    hash;    // hash    // not necessary to store, can be quickly computed from ibox

	HashMapField(){
		object  = NULL;
		n       = 0;
		bucket  = 0; // FIXME this is not necessary
	}

	inline void set( TYPE1* object_, ULONG bucket_ ){
		object  = object_;
		bucket  = bucket_;
		//hash    = hash_;
	}

};

//#define INFINITE_LOOP_DEBUG(ss)  if( (i==hash)||(i==(hash-1)) ){ printf( " INFINITE LOOP !!! %s \n ", ss );  return -i; };

#define INFINITE_LOOP_DEBUG(ss)

template <class TYPE >
class HashMap{
	public:
	UINT power;
	UINT mask;
	UINT capacity;
	UINT filled=0;

	int DEBUG_counter;

	HashMapField<TYPE>*  fields;

	void init( UINT power_ ){
		power      = power_;
		capacity   = 1<<power;
		mask       = capacity-1;
		//printf( "  power %i capacity %i mask %i \n ", power, capacity, mask );
		fields     = new HashMapField<TYPE>[ capacity ];
		//for( int i=0; i<capacity; i++ ){ fields[i].set( 0, 0 ); }   // initialization done by constructor (?)
	}

	inline void insert( TYPE* object, ULONG bucket, UINT hash, UINT i ){
		filled++;
		fields[ hash ].n++;
		fields[ i    ].set( object, bucket );
	}
	inline void remove( UINT i, UINT hash ){
		filled--;
		fields[ hash ].n--;
		fields[ i    ].object = NULL;
	};
	inline int findFreePlace( UINT hash ) const {
		UINT i    = hash;
		while( fields[i].object != NULL ){
			i=(i+1)&mask;
			INFINITE_LOOP_DEBUG( "findFreePlace" )
		}
		return i;
	}
	inline UINT insertNoTest( TYPE* object, ULONG bucket ){
		//if ( (filled<<1) > capacity ) resize( power+1 );  // FIXME : THIS COULD BE VERY UNSAFE
		UINT hash = mask&hashFunc( bucket );   // must be after resize
		UINT i    = findFreePlace( hash );
		insert( object, bucket, hash, i );
		return i;
	};

	void resize( UINT power_ ){
		UINT old_capacity = capacity;
		HashMapField<TYPE>*  old_fields = fields;
		init( power_ );
		HashMapField<TYPE>*  p_field    = old_fields;
		for (int i=0; i<old_capacity; i++){
			insertNoTest( p_field->object, p_field->bucket );
			p_field++;
		}
		delete old_fields;
	}

	inline int find( TYPE* object, ULONG bucket, UINT hash ) const {
		UINT n    = fields[ hash ].n;
		UINT i    = hash;
		//printf( " find : object %i bucket %i hash %i n %i \n", object, bucket, hash, n );
		while( n > 0 ){
			//printf( " : i %i n %i bucket %i hash(bucket) %i hash %i object %i \n", i, n, fields[i].bucket, mask&hashFunc( fields[ i ].bucket ), hash, fields[i].object );
			if( object == fields[ i ].object             ){ return i; }
			if( hash   == ( mask&hashFunc( fields[ i ].bucket ) ) ){ n--; } // FIXME : here we call hashFunc in loop, would storing hash improve performance ?
			i=(i+1)&mask;
			INFINITE_LOOP_DEBUG("find")
		}
		return -1;
	}
	inline int find( TYPE* object, ULONG bucket ) const {
		UINT hash = mask&hashFunc( bucket );
		return find( object, bucket, hash );
	}
	inline int insertIfNew( TYPE* object, ULONG bucket ){
		UINT hash = mask&hashFunc( bucket );
		int  i    = find( object, bucket, hash );
		if( i < 0 ){
			if ( (filled<<1) > capacity ) resize( power+1 );
			UINT i = findFreePlace( hash );
			insert( object, bucket, hash, i );
			return i;
		}
		return -1;
	};


	bool tryRemove( TYPE* object, ULONG bucket ){
		UINT hash = mask&hashFunc( bucket );
		int  i    = find( object, bucket, hash );
		if( i > 0 ){
			remove( i, hash);
			return true;
		};
		return false;
	};

	UINT getBucketIndexes( ULONG bucket, UINT * outi )  const {
		UINT hash = mask&hashFunc( bucket );
		UINT n    = fields[ hash ].n;
		UINT i    = hash;
		UINT j    = 0;
		while( n > 0 ){
			ULONG bucketi = fields[ i ].bucket;
			if( hash == ( mask&hashFunc( bucketi ) ) ){   // FIXME : here we call hashFunc in loop, would storing hash improve performance ?
				if( bucketi == bucket ){
					outi[ j ] = i;
					j++;
				}
				n--;
			}
			i=(i+1)&mask;
			INFINITE_LOOP_DEBUG("getBucketIndexes")
		}
		return j;
	};

	UINT getBucketObjects( ULONG bucket, TYPE** out ) const {
		UINT hash = mask&hashFunc( bucket );
		UINT n    = fields[ hash ].n;
		UINT i    = hash;
		UINT j    = 0;
		//printf( " getBucketObjects: bucket (%i,%i) hash %i n %i \n", (bucket&0xFFFF),((bucket>>16)&0xFFFF), hash, n );
		while( n > 0 ){
			ULONG bucketi = fields[ i ].bucket;
			if( hash == ( mask&hashFunc( bucketi ) ) ){   // FIXME : here we call hashFunc in loop, would storing hash improve performance ?
				if( bucketi == bucket ){
					out[ j ] = fields[ i ].object;
					j++;
				}
				n--;
			}
			i=(i+1)&mask;
			INFINITE_LOOP_DEBUG("getBucketObjects")
		}
		return j;
	};

};


#endif

