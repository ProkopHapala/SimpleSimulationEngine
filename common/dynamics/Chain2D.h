class Chain2D{
	public:
	int n;
	Vec2d  * poss;
	double * ls;


	void allocate( ){
		poss = new Vec2d [n];
		ls   = new double[n];
	}

	Chain2D( int n_ ){
		n = n_;
		allocate();
	}

	Chain2D( int n_, const Vec2d& p0, const Vec2d& d ){
		n = n_;
		allocate( );
		Vec2d p;
		p.set( p0 );
		double l = d.norm();
		for( int i=1; i<n; i++){
			p.add( d );
			poss[i].set( p );
			ls[i] = l;
		}
	}

	void move( const Vec2d& p0 ){
		Vec2d p;	
		p.set( p0 );
		for( int i=0; i<n; i++){
			Vec2d d,p_;
			p_.set( poss[i] );
			d.set_sub( p_, p );
			double r2 = d.norm2();
			double l = ls[i];
			if( r2 > l*l ){
				double r = sqrt(r2);
				p_.add_mul( d, l - r );
				poss[i-1].set( p_ );
			} 
			p.set( p_ );
		}
	}

	void draw( const Vec2d& p0  ){
		glBegin( GL_LINES );
			Vec2f p;
			convert(p0, p );
			for( int i=0; i<n; i++){
				Vec2f p_;
				convert( poss[i], p_ );
				glVertex3f( p .x, p .y, 0.0 );
				glVertex3f( p_.x, p_.y, 0.0 );
				p.set(p_);
			}
		glEnd();
	}

};
