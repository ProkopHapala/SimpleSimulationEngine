#ifndef HashMap64_h
#define HashMap64_h

/*
NOTES:

Size:
 - 64 bit boxId is 21 bit/axis = 2,097,152;   with 1km per box this is just 2 milion km per axis, which does not even cover Jupiter system
 - 32 bit per axis will allow 4 billion ticks ( 4 billion km per axis with 1km resolution ) that is approx distance to pluto
 - 64 bit per axis will allow 18.446744 billion kilometers (1.949 ly ) with 1 mm resolution  (light year 9.4607e+15 m)
 - 52 bit mantisa of double means 4.5e+15 ticks ( 0.5 ly with 1 m resolution ), distance to Pluto (4.43682e+9 km) with 1mm resolution
 - anyway trepasing >1000 boxes per raycast is super inefficient

*/


class HashMap64{ public:
	uint8_t  power    = 0;
	uint32_t capacity = 0;
	uint32_t mask     = 0;
	uint32_t filled   = 0;
	float    maxFillRatio = 0.6;
	//float    begFillRatio = 0.3;

	int DEBUG_counter;

	uint64_t EMPTY_I = 0xFFFFFFFFFFFFFFFFL;
	uint64_t EMPTY_P = 0xFFFFFFFFFFFFFFFFL;
	uint32_t EMPTY_H = 0xFFFFFFFF;

	uint64_t*  store;  // pointer or index of stored object
	uint64_t*  iboxs;  // unique box index
	uint32_t*  hashs;  // hash   // not necessary to store, can be quickly computed from ibox
	uint16_t*  hits;   // number of objects with this hash
	//uint8_t*  hits;   // number of objects with this hash

	//BucketId* iboxs;  // unique box index - we may need more than 64 bit number if we have large axis

	inline uint32_t hash( uint64_t ibox ){
		// see http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
		// https://stackoverflow.com/questions/6943493/hash-table-with-64-bit-values-as-key
		//return  ( 2166136261UL ^ i * 16777619 ) & mask;
		//return  ( 2166136261UL ^ (i * 16777619) ) & mask;
		//i = 1664525*(i+15445) ^ 1013904223; return (1664525*i ^ 1013904223) & mask;
		//return ((i >> power)^i) & mask;
		//return (i*2654435761 >> 16)&mask;   // Knuth's multiplicative method
        //printf( "hash ibox %li ", ibox );
		//return (ibox*2654435761 >> 16)&mask;   // Knuth's multiplicative method
		//return ( ( (ibox * 11400714819323198549UL) >> 32 ) + ibox )&mask;  // Knuth's multiplicative method  64bit
		//return ( ( (ibox * 11400714819323198549UL) >> 32 )  ^  ( ( (ibox>>32) * 11400714819323198549L) >> 32 )  )&mask;
		return ( ( ( ibox * 11400714819323198549UL )  ^  ( (ibox>>32) * 11400714819323198549UL ) ) >> 32 )&mask;
		//return ( ( ( ((uint32_t)ibox)*265443576) ^ ( (uint32_t)(ibox>>32)*265443576) ) >> 16 )&mask;
		//return ( ibox )&mask;  // Knuth's multiplicative method  64bit
		//i = ((i >> 16) ^ i) * 0x45d9f3b; i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
		//i = ((i >> 16) ^ i) * 0x3335b369; i = ((i >> 16) ^ i) * 0x3335b369; return  ((i >> 16) ^ i)&mask;
	}

	inline void set( int i, uint64_t p, uint64_t ibox, uint32_t h ){
		store[ i ] =  p;
		iboxs[ i ] =  ibox;
		hashs[ i ] =  h;
	}

    void clear(){
        filled     = 0;
        for (int i=0; i<capacity; i++){
			hits[i] =  0;
			set( i, EMPTY_P, EMPTY_I, EMPTY_H );
		}
	}

	void init( int power_ ){
		power      = power_;
		capacity   = 1<<power;
		mask       = capacity-1;
		//printf( "  power %i capacity %i mask %i \n ", power, capacity, mask );
		store      = new uint64_t[capacity];
		hits       = new uint16_t[capacity];
		hashs      = new uint32_t[capacity];
		iboxs      = new uint64_t[capacity];
        clear();
	}

	inline int getAllInBox( uint64_t ibox, uint32_t h, int* outi ){
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
			}
			i=(i+1)&mask;
			DEBUG_counter++;
		}
		return j;
	}
	inline int getAllInBox( uint64_t ibox, int* outi ){ uint32_t h = hash( ibox ); return getAllInBox( ibox, h, outi ); }

    inline int getAllInBoxOnce( uint64_t ibox, uint32_t h, int* outi, bool* isOld ){
		DEBUG_counter =0;
		int n = hits[h];
		int i = h;
		int j =0;
		//printf( "getAllInBoxOnce i: %i\n", i);
		while( n>0 ){
            //printf( "getAllInBoxOnce n: %i\n", n);
			if( hashs[i]==h ){
				if( iboxs[i]==ibox ){
                    uint32_t io = (uint32_t)store[i]; // Take Low bytes
                    //printf( "io: %i \n", io);
                    if( !isOld[io] ){
                        isOld[io] = true;
                        outi[j]   = i;
                        j++;
					}
				}
				n--;
			}
			i=(i+1)&mask;
			//printf( "getAllInBoxOnce i: %i\n", i);
			DEBUG_counter++; if(DEBUG_counter>capacity){  printf("getAllInBoxOnce: infinite loop \n"); exit(0); }
		}
		return j;
	}
	inline int getAllInBoxOnce( uint64_t ibox, int* outi, bool* isOld ){ uint32_t h = hash( ibox ); return getAllInBoxOnce( ibox, h, outi, isOld ); }


	inline int getAllInBox_noHash( uint64_t ibox, int* outi ){
		DEBUG_counter =0;
		uint32_t h = hash( ibox );
		int n = hits[h];
		int i = h;
		int j =0;
		while( n>0 ){
			int ib = iboxs[i];
			int ih = hash( ib ); // The difference is here
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

	inline void insert( uint64_t p, uint64_t ibox, uint32_t h, int i ){
		hits   [h] ++;
		set( i, p, ibox, h );
		filled++;
	}
	inline int insert( uint64_t p, uint64_t ibox, uint32_t h ){
        //printf( "insert( p %li ibox %li h %li \n", p, ibox, h );
        checkResize();
		int i = h;
		while( store[i] != EMPTY_P ) i=(i+1)&mask;
		insert( p, ibox, h, i );
		return i;
	};
	inline int insert( uint64_t p, uint64_t ibox ){ uint32_t h = hash( ibox ); return  insert( p, ibox, h );   }

	void remove( int i, uint32_t h ){
		filled--;
		hits[ h ]--;
		set( i, EMPTY_P, EMPTY_I, EMPTY_H );
	};
	void remove( int i ){ uint32_t h = hash( iboxs[i] ); remove( i, h ); };

	void resize( int power_ ){
		int old_capacity = capacity;
        uint64_t*  old_store = store;
        uint64_t*  old_iboxs = iboxs;
        //uint32_t*  old_hashs = hashs;
        if( old_capacity > 0 ){
            delete [] hits;
            delete [] hashs;
        }
		init( power_ );
		if( old_capacity > 0 ){
            for (int i=0; i<old_capacity; i++){
                if( old_store[i] != EMPTY_P ){
                    insert( old_store[i], old_iboxs[i] );
                }
            }
            delete [] old_store;
            //delete [] old_hashs;
            delete [] old_iboxs;
		}
		printf( "HashMap64:resized power %i capacity %i filled %i ratio %g \n", power, capacity, filled, filled/(float)capacity );
	}

	void checkResize(){
	   if( filled > capacity * maxFillRatio ){
	        //power++;
            while( ((1<<power)*maxFillRatio) < filled ) power++;
            //printf( "checkResize n %i n*f %i nMaxNew %i pow %i \n", filled, (int)(filled/begFillRatio), 1<<power, power );
	        resize( power );
	   };
	}

	int reserve( int n ){
        int nTarget = (int)(n/maxFillRatio)+1;
        //printf( "nTarget  %i \n", nTarget );
        if( nTarget > capacity ){
            while( (1<<power) < nTarget ){
                //printf( "power %i n %i \n", power, (1<<power) );
                power++;
                //if(power>16) break;
            }
            resize( power );
        }
        return capacity;
	}

	inline int nbytes(){ return sizeof(*this) + ( capacity*(  sizeof(uint64_t)*2 +  sizeof(uint32_t) + sizeof(uint16_t)  ) ); }

};

#endif
