

## Image

```glsl
#define T(ch)    texelFetch(ch,ivec2(U) % textureSize(ch,0),0).x
#define hash(p)  fract(sin(dot(p, vec2(11.9898, 78.233))) * 43758.5453) // iq suggestion, for Windows
// #define hash(p)  fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453) // classical verison

// see power spectrum here: https://www.shadertoy.com/view/7dyfzW
#if 1
float B(vec2 U) {
    float v = 0.;
    for (int k=0; k<9; k++)
        v += hash( U + vec2(k%3-1,k/3-1) ); 
  //return       1.125*hash(U)- v/8.  + .5; // some overbound, closer contrast
    return .9 *( 1.125*hash(U)- v/8.) + .5; // 
  //return .75*( 1.125*hash(U)- v/8.) + .5; // very slight overbound
  //return .65*( 1.125*hash(U)- v/8.) + .5; // dimmed, but histo just fit without clamp. flat up to .5 +- .23
}
#else
float B(vec2 U) {                           // 5-tap version 
    float v =  hash( U + vec2(-1, 0) )
             + hash( U + vec2( 1, 0) )
             + hash( U + vec2( 0, 1) )
             + hash( U + vec2( 0,-1) ); 
    return  hash(U) - v/4.  + .5;
}
#endif

void mainImage(out vec4 O, vec2 U) {
    
    float T = iTime, t = iMouse.z<=0. ? 0. : float(77*iFrame);
    vec3  R = iResolution,
          D = vec3(.3*(U+U-R.xy)/R.y, -1),  // ray direction
          p = 30./R, q;                     // marching point along ray 
    O-=O;
    int x = int(4.*U.x/R.x);
    U = mod(U-t,R.xy);
    if (D.y<0.) U = floor(U/2.);
    O +=   x == 0 ? T(iChannel0)            // white texture noise
         : x == 1 ? T(iChannel1)            // blue texture noise
         : x == 2 ? hash(U)                 // white  procedural noise
         :          B(U);                   // blue procedural noise
    
  //O.x>1. ? O-=O,O.x++ : O.x<0. ? O-=O,O.y++ : t ;   // verify not out of bounds
}
```
