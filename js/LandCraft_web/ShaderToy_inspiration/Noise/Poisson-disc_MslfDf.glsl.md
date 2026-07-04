
## Buffer A

```glsl
// Poisson disk : repealing particles.

#define rnd(U)  ( 2. * fract(4567.89*sin(4567.89*U*mat2(1,-13.17,377.1,-78.73))) - 1. )
 
#define T(i,j)  texture(iChannel0, fract( ( U + vec2(i,j) ) / iResolution.xy ) ).xy

void mainImage( out vec4 O, vec2 U )
{
    if (iFrame==0) {
        O.xy = .2*rnd(U);                   // particles location = vec2(i,j) + stored perturb
        return;                             // start with jittered grid (rough approx of Poisson disc)
    }
    
    vec2 U0 = T(0,0), D, F = vec2(0);
    
    for (int j=-4; j<=4; j++)               // look in the neighborhood
        for (int i=-4; i<=4; i++)           // with [-3,3]^2, a leak occurs at ~100"+resize
            if ( vec2(i,j) != vec2(0) ) {
                D = vec2(i,j)+T(i,j) - U0;  // distance vector to particle (i,j)
                float l = length(D);
                F += D / l * max(2.-l,0.);  // simulates a spring (only repealing, otherwise clamped)
            }
    
    O.xy = U0 - .1* F;                      // displace particle proportionaly to force
}
```

## Image

```glsl
// Generate Poisson-disc distribution by simulated repulsive springs between particles.
// see spectrum (and tilable version) in https://www.shadertoy.com/view/MssfDf

#define zoom 10.

#define P(i,j)  ( floor(U/zoom)+vec2(i,j)+ texelFetch(iChannel0, ivec2(U/zoom)+ivec2(i,j), 0).xy ) 

//#define N 30
//#define P(i,j)  ( vec2(i,j) + texelFetch(iChannel0, ivec2(i,j), 0).xy )  // for 0..N loops


/**/
void mainImage( out vec4 O, vec2 U )
{   
    O -= O; 
 
    for (int j=-3; j<=3; j++)
        for (int i=-3; i<=3; i++) 
            O += smoothstep(2.,0., length( P(i,j) *zoom - U ) );
}
/**/




/**  // Voronoi version

void mainImage( out vec4 O, vec2 U )
{   
    float v, l = 1e9,_l;
    
    for (int j=-3; j<=3; j++)
        for (int i=-3; i<=3; i++) {
            v = length( P(i,j) *zoom - U );
            if (v<l) _l=l, l=v;
        }
    
    O = .1* vec4(_l-l);  
    //O = .1*vec4(l,0,_l-l,0);
}
/**/

```
