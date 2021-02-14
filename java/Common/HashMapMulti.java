package Common;


///
/// Implementation of Open-addresing hasm map storing multiple elements per cell
/// https://www.geeksforgeeks.org/hashing-set-3-open-addressing/


public class HashMapMulti{
    public int power;
    public int mask;
    public int capacity;
    public int filled=0;
    public int DEBUG_counter;

    public Object[] store;
    public int[]    hits;  // number of objects with this hash
    public long[]   iboxs; // unique box index (before hashing, independent of hash-map size)
    //public int[]  hashs; // hash   // not necessary to store, can be quickly computed from ibox

    public final int hash( long i ){
        // see http://stackoverflow.com/questions/664014/what-integer-hash-function-are-good-that-accepts-an-integer-hash-key
        //return  ( 2166136261UL ^ i * 16777619 ) & mask;
        //return  ( 2166136261UL ^ (i * 16777619) ) & mask;
        //i = 1664525*(i+15445) ^ 1013904223; return (1664525*i ^ 1013904223) & mask;
        //return ((i >> power)^i) & mask;
        i=((int)((i*2654435761L) >> 16));
        return ((int)((i*2654435761L) >> 16))&mask;   // Knuth's multiplicative method
        //i = ((i >> 16) ^ i) * 0x45d9f3b; i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
        //i = ((i >> 16) ^ i) * 0x45d9f3b; return  ((i >> 16) ^ i)&mask;
        //i = ((i >> 16) ^ i) * 0x3335b369; i = ((i >> 16) ^ i) * 0x3335b369; return  ((i >> 16) ^ i)&mask;
    }

    public final void set( int i, Object p, long ibox ){
        store  [ i ] =  p;
        iboxs  [ i ] =  ibox;
        //hashs  [ i ] =  h;
    }
    
    public final void checkResize(){
        //System.out.println( "  filled "+(filled<<1)+" capacity "+capacity  );
        if ( (filled<<1) > capacity ) resize( power+1 );   // fill-factor = 0.5
    }
    
    public final void put( int i, Object p, long ibox, int h ){
        checkResize();
        //System.out.println(  "insert ibox: "+ibox+" h: "+h+" i: "+i+" n[h] "+hits[h] );
        hits[h]++;
        set ( i, p, ibox );
        filled++;
    }
    
   public final int insert( Object p, long ibox ){
        int h = hash( ibox );
        int i = h;
        while( store[i] != null  ){ i=(i+1)&mask; }
        put( i, p, ibox, h );
        return i;
    };

    public final void remove( int i ){
        filled--;
        int h = hash( iboxs[i] );
        hits[ h ]--;
        set( i, null, -1 );
    };
    
    public final int getAllInBox( long ibox, Object[] outi ){
        //DEBUG_counter =0;
        int h = hash( ibox );
        int n = hits[h];
        //System.out.println(  "getAllInBox n: "+n );
        int i=h;
        int j =0;
        while( n>0 ){
            long ib = iboxs[i];
            int  ih = hash( ib );
            if(  ih==h ){
                if( ib==ibox ){
                    //outi[j] = i;
                    outi[j]=store[i];
                    j++;
                }
                n--;
            }
            i=(i+1)&mask;
            //DEBUG_counter++;
        }
        //System.out.println( "hits[h] "+hits[h]+" : "+j  );
        return j;
    }

    public final void resize( int power_ ){
        //int old_filled      = filled;
        int old_capacity  = capacity;
        Object[]  old_store = store;
        //int[]   old_hashs = hashs;
        long[]    old_iboxs = iboxs;
        init( power_ );
        //if(capacity>100)System.exit(0);
        for (int i=0; i<old_capacity; i++){
            if( old_store[i] !=null ){
                insert( old_store[i], old_iboxs[i] );
            }
        }
    }
    
   public final void init( int power_ ){
        power      = power_;
        capacity   = 1<<power;
        mask       = capacity-1;
        //System.out.println( "  power "+power+" capacity "+capacity+" mask "+ mask  );
        //System.out.println( String.format( "  power %i capacity %i mask %i ", power, capacity, mask ) );
        store      = new Object[capacity];
        hits       = new int   [capacity];
        //hashs    = new int   [capacity];
        iboxs      = new long  [capacity];
        filled=0;
        for (int i=0; i<capacity; i++){
            hits  [i] =  0    ;
            set( i, null, -1 );
        }
    }
   
   
    public HashMapMulti( int power_ ){
       init( power_ );
   }
    
}