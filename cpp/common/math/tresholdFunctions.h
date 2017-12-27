
#ifndef  tresholdFunctions_h
#define  tresholdFunctions_h

namespace Treshold{

	// ==== polynom order 0

	inline double p0i( double x ){
		if  ( x < 0 ){ return 0; }
		else         { return x; }
	};

	inline double p0( double x ){
		if  ( x < 0 ){ return 0; }
		else         { return 1; }
	};

	// ==== polynom order 1

	inline double p1i( double x ){
		if      ( x < 0 ){ return 0;        }
		else if ( x > 1 ){ return x-0.5d;   }
		else             { return 0.5d*x*x; }
	};

	inline double p1( double x ){
		if      ( x < 0 ){ return 0; }
		else if ( x > 1 ){ return 1; }
		else             { return x; }
	};

	inline double p1d( double x ){
		if      ( x < 0 ){ return 0; }
		else if ( x > 1 ){ return 0; }
		else             { return 1; }

	};

	// ==== polynom order 2 // two parabolas

	inline double p2i( double x ){
		if      ( x < 0   ){ return 0; }
		else if ( x > 1.0 ){ return x-0.5d; }
		else if ( x < 0.5 ){ return    0.66666666666*x*x*x; }
        else               { return ((-0.66666666666*x + 2)*x - 1 )*x + 0.16666666666; }
		//else             { return -0.66666666666*x*x*x + 2*x*x - x + 0.16666666666; }
	};

	inline double p2( double x ){
		if      ( x < 0   ){ return 0;    }
		else if ( x > 1.0 ){ return 1;    }
		else if ( x < 0.5 ){ return 2*x*x; }
		else               { return (4-2*x)*x-1; }
        //else               { double x_ = x-1; return 1-2*x_*x_; }
		//else             { return 1-2*x*x+4*x-2; }
	};

	inline double p2d( double x ){
		if      ( x < 0   ){ return 0;        }
		else if ( x > 1.0 ){ return 0;        }
		else if ( x < 0.5 ){ return 4*x;      }
		else               { return 4*(1-x);  }
	};

	// ==== polynom order 3 // cubic see https://en.wikipedia.org/wiki/Smoothstep

	inline double p3i( double x ){
		if      ( x < 0   ){ return 0;                       }
		else if ( x > 1.0 ){ return x - 0.5d;                   }
		else               { return x*x*x*( 1 - 0.5d * x );  }
	};

	inline double p3( double x ){
		if      ( x < 0   ){ return 0;           }
		else if ( x > 1.0 ){ return 1;           }
		else               { return x*x*(3-2*x); }
	};

	inline double p3d( double x ){
		if      ( x < 0   ){ return 0;          }
		else if ( x > 1.0 ){ return 0;          }
		else               { return 6*x*(1-x);  }
	};


    // ==== polynom order 4 // NOT GOOD
    /*
	inline double p4i( double x ){
		if      ( x < 0   ){ return 0; }
		else if ( x > 1.0 ){ return x-0.5d; }
		else if ( x < 0.5 ){ return    0.66666666666*x*x*x; }
        else               { return ((-0.66666666666*x + 2)*x - 1 )*x + 0.16666666666; }
		//else             { return -0.66666666666*x*x*x + 2*x*x - x + 0.16666666666; }
	};

	inline double p4( double x ){
		if      ( x < 0   ){ return 0;    }
		else if ( x > 1.0 ){ return 1;    }
		else if ( x < 0.5 ){ double x2 =x*x;                     return     8*x2*x2; }
		else               { double x_ = x-1; double x2 = x_*x_; return 1 - 8*x2*x2; }
        //else               { double x_ = x-1; return 1-2*x_*x_; }
		//else             { return 1-2*x*x+4*x-2; }
	};

	inline double p4d( double x ){
		if      ( x < 0   ){ return 0;        }
		else if ( x > 1.0 ){ return 0;        }
		else if ( x < 0.5 ){ return  32*x*x*x;      }
		else               { return -32*x*x*x;  }
	};
	*/


	// ==== polynom order 5 // Perlin see https://en.wikipedia.org/wiki/Smoothstep

	inline double p5i( double x ){
		if      ( x < 0   ){ return 0;                                  }
		else if ( x > 1.0 ){ return x-0.5d;                             }
		else               { return x*x*x*x*( 10.0d/4 - x*( 3 - x ) );  }
	};

	inline double p5( double x ){
		if      ( x < 0   ){ return 0;                            }
		else if ( x > 1.0 ){ return 1;                            }
		else               { return x*x*x*(10 - x*( 15 - 6*x ) ); }
		//else               { double x2=x*x; return x*x2*( x2*6 - 15*x +10 ); }
	};

	inline double p5d( double x ){
		if      ( x < 0   ){ return 0;          }
		else if ( x > 1.0 ){ return 0;          }
		else               { return x*x*( 30 - x * ( 60 - 30*x ) );  }
	};


	// ==== 1/x

    inline double r1i( double x ){
		if      ( x < 0   ){ return - x - log(1-x);  }
		else               { return   x - log(1+x);  }
	};


    inline double r1( double x ){
		if      ( x < 0   ){ return x/(1-x);  }
		else               { return x/(1+x);  }
	};

    inline double r1d( double x ){
		if      ( x < 0   ){ double x_ = 1-x; return 1/(x_*x_); }
		else               { double x_ = 1+x; return 1/(x_*x_); }
	};


    // ==== 1/x2

    inline double r2i( double x ){
		if      ( x < 0   ){ double x_ = 1-x; return -x + 0.5d/(x_*x_) - 0.5d;   }
		else               { double x_ = 1+x; return  x + 0.5d/(x_*x_) - 0.5d;  }
	};


    inline double r2( double x ){
		if      ( x < 0   ){ double x_ = 1-x; return -1 + 1/(x_*x_*x_);  }
		else               { double x_ = 1+x; return  1 - 1/(x_*x_*x_);  }
	};

    inline double r2d( double x ){
		if      ( x < 0   ){ double x_ = 1-x; double x2 = x_*x_; return 3/(x2*x2); }
		else               { double x_ = 1+x; double x2 = x_*x_; return 3/(x2*x2); }
	};


    // ==== 1/sqrt(1+x2)

    inline double root2i( double x ){
        return sqrt( 1 + x*x );
	};


    inline double root2( double x ){
        return x/sqrt( 1 + x*x );
	};

    inline double root2d( double x ){
        double a = sqrt(1 + x*x);
        return 1/(a*a*a);
	};

}

#endif
