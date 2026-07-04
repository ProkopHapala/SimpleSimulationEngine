

## Image

```glsl
// https://en.wikipedia.org/wiki/Lanczos_resampling

const float a = 3., n = a-.5, N = 2., PI = 3.1415927;

//#define L(i)    ( i==0. ? 1. : a* sin(PI*i)*sin(PI*i/a) / ( PI*PI*i*i ) )
  #define L(i) float[]( 1., .60793, 0., -.13509, 0., .024317 )[abs(int(i+i))] // for a=3,N=2

#define T(i,j)  ( 2.* texture( iChannel0, ( u +N*vec2(i,j) ) / 256. ) - 1. )

#define blue() (  \
          (  H(-1,-1) +     H(0,-1)  + H(1,-1)  \
           + H(-1, 0) - 8.* H(0, 0)  + H(1, 0)  \
           + H(-1, 1) +     H(0, 1)  + H(1, 1)  \
          ) *.5/9. *2.1 +.5 )
#define H(i,j) texture( iChannel0, ( u + vec2(i,j) ) /256. )


void mainImage( out vec4 O, vec2 u )
{
    vec2 R = iResolution.xy,
         U = u / R;
    bool b = iMouse.z<=0.;
    if ( b && int(u) % int(R/3.) == 0 ) { O = vec4(1,0,0,1); return; } // red separator
    int  i = int( 3.*( b ? u : iMouse.xy ) / R );  // select noise via zone or mouse
    
    if ( b && U.y < .5 ) u = floor(u/2.)+.5;       // bottom: zoom
    
    if ( i==0 )                                    // left: cheap blue with 9 points FD
        O = blue();
    
    else if (i==1) {                               // middle: blue with 11x11  Lanczos filtering
        O *= 0.;                              
        for( float y=-n; y<=n; y+= 1./N )
            for( float x=-n; x<=n; x+= 1./N )
                O +=   L(x)*L(y) * T(x,y);

        O = .5 + .65*( T(0,0) - O/(N*N) );         // keep only high freqs
                                                   // testing some histogram corrections
     // O = .5 + .7*tanh((O*2.-1.)/(1.76*.8));     // would be valid for gaussian distrib
     // O = (O-.5)*.6+.5;                          // fit to [0,1] interval → histo = /‾\
        // O = 2.*O-1.; O = .5+.5*(O + sin(1.57*O))/2.;
        // O = (O + O*O*(3.-2.*O))/2.;
        // float s = .242, x = abs(O.x-.5), k = 3.14/(.5-s); O = vec4(.5+.5*sign(O.x-.5)*.517*( x<s ? s*x : s*s+ .5*x+.5/k*sin(k*(x-s))));
    }
    else                                           // right: ref blue noise
         O = texelFetch(iChannel1, ivec2(u)%1024,0);

    if ( b && abs( U.y-.5) < .25 ) O = step(.8,O); // centerline: thresholding
    O = O.xxxx;                                    // grey
}
```
