

## Image

```glsl
// various variants             https://www.shadertoy.com/view/tllcR2
// LOD0-LOD1 on 9 points + FFT  https://www.shadertoy.com/view/7dyfzW
   
#define H(p)  fract(sin(mod(dot(p, vec2(12.9898, 78.233)),6.283)) * 43758.5453)

#define blue(p) (  \
          (  H(p+vec2(-1,-1)) + H(p+vec2(0,-1)) + H(p+vec2(1,-1))  \
           + H(p+vec2(-1, 0)) - 8.* H( p )      + H(p+vec2(1, 0))  \
           + H(p+vec2(-1, 1)) + H(p+vec2(0, 1)) + H(p+vec2(1, 1))  \
          ) *.5/9. *2.1 +.5 )
   
#define T(l) texelFetch(iChannel0, (ivec2(u)%256) >> l, l )

void mainImage( out vec4 O, vec2 u )
{
    u -= .5;
    if ( int(u) % int(iResolution/4.) == 0 ) { O = vec4(1,0,0,1); return; } // red separator
    vec2 U = u/iResolution.xy;
    int i = int(4.*U.x); 
    if (U.y<.5) u = floor(u/2.);          // bottom: zoom
    
    O =  i==0 ? vec4(H(u))                // left to right: white, cheapblue1, 2, ref blue
       : i==1 ?.5 + (T(0) -  T(1))
       : i==2 ? vec4(blue(u))
       :        texelFetch(iChannel1, ivec2(u)%1024,0);
    
    if (abs(U.y-.5)<.25) O = step(.9,O);  // center line: thresholding
    
    O = O.rrrr;                           // grey
}
```
