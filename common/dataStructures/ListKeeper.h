

class ListKeeper{
	public:
	int n;
	int * data;

	ListKeeper( int n_ ){
		n = n_;
		data  = new int[n];
		for( int i=0;i<n; i++){ data[i]=0; }
	};

	void grow( int n_ ){
		int * data_  = new int[n_];
		for( int i=0;i<n; i++){ data_[i]=data[i]; }
		data = data_;
		n = n_;
	};

	int getFree(){
		int freelits = -1;
		for( int i=0; i<n; i++ ){ if( data[i]==0 ) freelits=i; break; }
		if( freelits == -1 ){
			freelits = n;
			grow( (int)( 1 + n*1.61803398875 ) );
		}
		return freelits;
	};

	int countEmpty( ){
		int empty = 0;
		for( int i=0;i<n; i++){ if( data[i]==0 ) empty++; }
		return empty;
	};

	void shrink( int n_ ){
		int * data_  = new int[n_];
		int j = 0;
		for( int i=0;i<n; i++){ 
			if( data[i] > 0 ){
				data_[i]=data[i];
				j++; 
			}			
		}
		data = data_;
		n = n_;
	}

	bool tryShrink( float triggerFill , float targetFill ){
		int empty = countEmpty();
		if( (empty/(float)n) < triggerFill ){
			shrink( (int)(n*targetFill+1) );
			return true;
		}
		return false; 
	}

};
