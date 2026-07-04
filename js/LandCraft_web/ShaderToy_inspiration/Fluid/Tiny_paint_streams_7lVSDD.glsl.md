
## Common

```glsl
#define R iResolution.xy
#define A texture(iChannel0,vec2(u+i)/R)
#define M void mainImage(out vec4 r, vec2 u) { vec2 x, i = u-u, f = i; r -= r
```

## Buffer A

```glsl
M; float G, m = A.w;
    // iterate over neighbouring cells
    // A.xy = velocity
    // A.z = pigment
    // A.w = density
    for(int k = 81; k-->0;)
        i = vec2(k%9,k/9)-4.,
        x = i + A.xy, // move particle according to velocity

        // Gaussian diffusion, sigma=sqrt(.5)
        G = A.w / exp(dot(x,x)) / 3.142,

        // advection
        r += vec4(A.xyz, 1) * G,
    
        // pressure forces (smoothed particle hydrodynamics)
        f -= ( m*m-m       // pressure at current position
                           // density * (density - reference fluid density)
             + A.w*A.w-A.w // pressure at neighbour position
             ) * G * x;    // gradient of smoothing kernel

    if(r.w > 0.) // not vacuum
        r.xyz /= r.w, // convert momentum to velocity, normalise pigment
        r.xy += clamp(f / r.w, -.1, .1); // acceleration
    
    // gravity
    r.y -= .005;

    // boundary
    if(u.y < 9.) r.y += 1.;
    
    // streams
    r = length(u - R * vec2(.2, .9)) < 9. ?
            vec4(sin(iTime) + 2., -1, 0, 1) :
        length(u - R * vec2(.8, .9)) < 9. ?
            vec4(sin(iTime) - 2., -1, 1, 1) :
        iFrame < 2 && u.y > .8*R.y ? // init
            vec4(0,0, u.x < R.x/2., 1) :
        r;
}
```


## Image

```glsl
M - 1. + A.w * vec4(0, A.z, 1. - A.z, 0); }
```
