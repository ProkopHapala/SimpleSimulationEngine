
## Common

```glsl
#define M void mainImage(out vec4 r, vec2 u)
#define A(i) texelFetch(iChannel0,ivec2(i+u),0)
```

## Buffer A

```glsl
M {
    r *= 0.;
    for(vec2 i = vec2(-7); ++i.x < 7.;) for(i.y = -7.; ++i.y < 7.;) {
        vec2 v = A(i).xy;                          // neighbour velocity
        r += A(i).z                                // neighbour mass
                * exp(-dot(v+i,v+i)) / 3.14        // normalised Gaussian
                * vec4(mix(v+v+i, v, A(i).z),1,1); // velocity contribution
    }
    r.xy /= r.z + 1e-6;
    if(iFrame % 500 == 1) {
        vec2 m = 4.*u/iResolution.xy-2.;
        r += vec4(m,1,0) * exp(-dot(m,m));
    }
}
```

## Image

```glsl
M {
    r = 1.-A().zzzz;
}
```
