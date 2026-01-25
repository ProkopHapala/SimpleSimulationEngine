

## Image

```glsl
// --- bluenoise:
// complicated ways:  "Stippling and Blue Noise" http://www.joesfer.com/?p=108 
// on shadertoy: https://www.shadertoy.com/view/4dK3WR
//               https://www.shadertoy.com/view/Md3GWf
//        my tests: https://www.shadertoy.com/view/ldyXDd

// --- perception of symmetries: 
// see https://www.shadertoy.com/view/XddSWn#

#define T(l) textureLod(iChannel0,U/256.,l)

float mask(vec2 p) { // see https://www.shadertoy.com/view/ldyXDd
#define DMUL  8.12235325   
#define SIZE  5.5
    vec2 U = floor(p/SIZE)*SIZE;
    p += ( T(0.).xy - .5 ) *DMUL;
    return fract( p.x*1.705 + p.y*.5375 ); 
}

void mainImage( out vec4 O, vec2 U )
{  
    float x = U.x/iResolution.x,
          y = U.y/iResolution.y;
    //U = floor(U/4.)+.5;
    
    if ( y < .33 ) {               // --- bottom: symmetries
          U = mod(U-=128.,255.);   // noise periodicity ( true 256 woul double lines )
        U.x = min(U.x,256.-U.x);   // noise symmetry in x
        U.y = min(U.y,256.-U.y);   // noise symmetry in y
    }
    
    O  =  x<.25 ? T(0.)
        : x<.26 ? O-O
        : x<.50 ? .5 + .5 * (T(0.)-.5) / T(1.)
        : x<.51 ? O-O 
        : x<.75 ? T(0.)-T(1.)+.5  //	O =( T(0.)-(T(1.)-.25) ) ;
        : x<.76 ? O-O
        :         vec4(mask(U));    // float(mask(U)>.8);

    y =  abs(y-.5)-.17;
    O = vec4( y>0. || O.x==0. 
                ? y < .01 ? 0.                // horizontal lines
                : O.x                         // top: homogeneous
                : step ( O.x, fract(4.*x)) ); // middle: ramp

}
```
