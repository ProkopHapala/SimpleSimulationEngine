
namespace Treshold{

	// ==== polynom order 0

	inline double p0d0( double x ){
		if  ( x < 0 ){ return 0; }
		else         { return x; }
	};

	inline double p0d1( double x ){
		if  ( x < 0 ){ return 0; }
		else         { return 1; }
	};

	// ==== polynom order 1

	inline double p2d0( double x ){
		if      ( x < 0 ){ return 0; }
		else if ( x > 1 ){ return 0.5d*x*x; }
		else             { return x+0.5d; }
	};

	inline double p1d1( double x ){
		if      ( x < 0 ){ return 0; }
		else if ( x > 1 ){ return 1; }
		else             { return x; }
	};

	inline double p1d2( double x ){
		if      ( x < 0 ){ return 0; }
		else if ( x > 1 ){ return 0; }
		else             { return 1; }

	};

	// ==== order 1 // two parabolas 

	inline double p2d0( double x ){
		if      ( x < 0   ){ return 0; }
		else if ( x > 1.0 ){ return x }
		else if ( x < 0.5 ){ return 0.66666666666*x*x*x;; }
		else               { return 1; }
	};

	inline double p2d1( double x ){
		if      ( x < 0   ){ return 0;    }
		else if ( x > 1.0 ){ return 1;    }
		else if ( x < 0.5 ){ return 2*x*x }
		else               { double x_ = x-1; return 1-2*x_*x_; }
	};

	inline double p2d2( double x ){
		if      ( x < 0   ){ return 0;        }
		else if ( x > 1.0 ){ return 1;        }
		else if ( x < 0.5 ){ return 4*x;      }
		else               { return 4*(1-x);  }
	};

	// ==== order 3 // cubic see https://en.wikipedia.org/wiki/Smoothstep

	inline double p3d0( double x ){
		if      ( x < 0   ){ return 0;                       }
		else if ( x > 1.0 ){ return 0.5+x;                   }
		else               { return x*x*x*( 1 - 0.5d * x );  }
	};

	inline double p3d1( double x ){
		if      ( x < 0   ){ return 0;           }
		else if ( x > 1.0 ){ return 1;           }
		else               { return x*x*(3-2*x); }
	};

	inline double p3d2( double x ){
		if      ( x < 0   ){ return 0;          }
		else if ( x > 1.0 ){ return 0;          }
		else               { return 6*x*(1-x);  }
	};


	// ==== Perlin order 5 // cubic see https://en.wikipedia.org/wiki/Smoothstep

	inline double p5d0( double x ){
		if      ( x < 0   ){ return 0;                       }
		else if ( x > 1.0 ){ return 0.5+x;                   }
		else               { return x*x*x*( 1 - 0.5d * x );  }
	};

	inline double p5d1( double x ){
		if      ( x < 0   ){ return 0;           }
		else if ( x > 1.0 ){ return 1;           }
		//else             { return x*x*x*(10 + x*( x*6 - 15 ) ); }
		else               { double x2=x*x; return x*x2*( x2*6 - 15*x +10 ) ); }
	};

	inline double p5d2( double x ){
		if      ( x < 0   ){ return 0;          }
		else if ( x > 1.0 ){ return 0;          }
		else               { return 6*x*(1-x);  }
	};

}
