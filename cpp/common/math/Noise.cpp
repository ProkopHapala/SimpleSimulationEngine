
#include "Noise.h" // THE HEADER

/*

  This is my atempt to produce the simplest and the fastest noise function in 2D
  it is composed of three parts:
   1) map given (x,y) coordinate to som trinagle (simplex) from regular simplex tiling
      we obtain index of triangle vertexes itx,ity and double coordinate inside the triangle (dx,dy)
   2) compute square of distance (r2) from the vertexes of the trinagle.
   3) Evaluate weights of vertexes as decaying spline. We don't want use sqrt() function because it is quite slow.
      function splineR2 is able to produce smooth function just from r2 without evaluation of r
   4) Evaluate pseudo-random hash from indexes of vertexes. We use multiply and xor. also other function and "magic numbers" can be used
   6) interpolate between the random values from the vertexes. Thanks to equally spaced vertexes we can approximate the interpolation by
      simple sum of weighted basisfuctions
*/

void Noise::simplexNoise2D( const Vec2d& pos, Vec2d& dpos ){
    double yf =  pos.y;
    double xf = (pos.x*0.86602540378 - 0.5*pos.y);
    int ity   = int( yf +10000.0 ) - 10000; // fast floor
    int itx   = int( xf +10000.0 ) - 10000;
    int seet  = (itx<<16)+ity;
    int seet2 = (itx<<16)+ity + 545645;
    double dx = itx*1.15470053839+ity*0.57735026919 - pos.x;
    double dy = ity                                 - pos.y;
    double w2 = splineR2( getR2( dx + 1.15470053839 , dy      ) );
    double w3 = splineR2( getR2( dx + 0.57735026919 , dy +1.0 ) );
    double cx = w2*rand_hash2( seet  + 0x10000 ) + w3*rand_hash2( seet  + 1 );
    double cy = w2*rand_hash2( seet2 + 0x10000 ) + w3*rand_hash2( seet2 + 1 );
    double w1;
    if( (yf + xf-(itx+ity) )<1.0 ){
        w1  = splineR2( getR2( dx , dy  ) );
        cx += w1 * rand_hash2( seet  );
        cy += w1 * rand_hash2( seet2 );
    }else{
        w1  = splineR2( getR2( dx + 1.15470053839+0.57735026919, dy + 1.0  ) );
        cx += w1 * rand_hash2( seet  + 0x10001 );
        cy += w1 * rand_hash2( seet2 + 0x10001 );
    }
    dpos.x = (cx/4294967296.0);
    dpos.y = (cy/4294967296.0);
    //printf( " (%f,%f) (%f,%f) \n", pos.x, pos.y, cx, cy );
}


void Noise::warpNoise3R( const Vec2d& pos0, const Vec2d& rot, double fdown, double strenght, int n, Vec2d& dpos_ ){
    Vec2d pos,dpos;
    pos  .set(pos0);
    dpos .set(0.0);
    for (int i=0; i<n; i++){
        //float xx = x*f+ sc*(noise_dx*ca + noise_dy*sa);
        //float yy = y*f+ sc*(noise_dy*ca - noise_dx*sa);
        dpos.mul_cmplx( rot );
        pos.mul( fdown );
        pos.add_mul( dpos, strenght );
        simplexNoise2D( pos, dpos );
    };
    dpos_.set(dpos);
    //dpos_.set_sub( pos, pos0 );
};
