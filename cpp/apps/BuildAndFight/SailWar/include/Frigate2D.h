
class Frigate2D : public Yacht2D {
	public:
	
	int nguns;
	Gun ** left_guns, right_guns; 

	void fire(){

	}

	virtual void draw( ){ 
		keel  .draw  ( *this );
		rudder.draw( *this );
		mast  .draw  ( *this );
	}

	bool loadFromFile( char const* filename ){
		//using namespace std;
		printf(" filename: %s \n", filename );
		FILE * pFile;
		pFile = fopen ( filename, "r" );
		const int nbuf = 1000; 
		char line [ nbuf ];
		fgets( line, nbuf, pFile );   keel  .fromString( line );  printf( "%s \n", keel  .toString( ) );
		fgets( line, nbuf, pFile );   rudder.fromString( line );  printf( "%s \n", rudder.toString( ) );
		fgets( line, nbuf, pFile );   mast  .fromString( line );  printf( "%s \n", mast  .toString( ) );
		fclose( pFile );
  		return false;
	}

};





