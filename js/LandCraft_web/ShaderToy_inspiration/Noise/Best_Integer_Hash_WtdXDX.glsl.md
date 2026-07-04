
## Image

```glsl
// --- from Chris Wellons https://nullprogram.com/blog/2018/07/31/
// Note that it might not be costlier than the infamous fract(big*sin(big*x)) ;-) 

        // --- choose one:
//#define hashi(x)   lowbias32(x)
  #define hashi(x)   triple32(x) 

  #define hash(x)  ( float( hashi(x) ) / float( 0xffffffffU ) )

//bias: 0.17353355999581582 ( very probably the best of its kind )
uint lowbias32(uint x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

// bias: 0.020888578919738908 = minimal theoretic limit
uint triple32(uint x)
{
    x ^= x >> 17;
    x *= 0xed5ad4bbU;
    x ^= x >> 11;
    x *= 0xac4c1b51U;
    x ^= x >> 15;
    x *= 0x31848babU;
    x ^= x >> 14;
    return x;
}

void mainImage( out vec4 O, vec2 U )
{
    uvec2 V = uvec2(U);
    float h = hash( V.x + hashi(V.y) ); // clean 2D hash
  //float h = hash( V.x + (V.y<<16) );  // 2D hash (should be ok too )
    O = vec4( h );
  //O = vec4( pow( h, 1./2.2) );        // sRGB conversion
}
```