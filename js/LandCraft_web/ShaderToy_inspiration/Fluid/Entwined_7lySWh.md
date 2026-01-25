
## Common

```glsl
/* Texture Stencil Library https://www.shadertoy.com/view/ssBczm

The MIT License

Copyright (c) 2022 David A Roberts <https://davidar.io/>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

vec4 texStencil(sampler2D ch, vec2 uv, float coeff[9]) {
    vec2 texel = 1. / vec2(textureSize(ch, 0));
    const vec2 stencilOffset[9] = vec2[](
        vec2(-1, 1), vec2( 0, 1), vec2( 1, 1),
        vec2(-1, 0), vec2( 0, 0), vec2( 1, 0),
        vec2(-1,-1), vec2( 0,-1), vec2( 1,-1)
    );
    vec4 r = vec4(0);
    for (int i = 0; i < 9; i++)
        r += coeff[i] * texture(ch, uv + texel * stencilOffset[i]);
    return r;
}

// Gaussian/binomial blur
// https://bartwronski.com/2021/10/31/practical-gaussian-filter-binomial-filter-and-small-sigma-gaussians/
vec4 texBlur(sampler2D ch, vec2 uv) {
    return texStencil(ch, uv, float[](
        .0625, .125, .0625,
        .125,  .25,  .125,
        .0625, .125, .0625
    ));
}

// Laplacian, optimal 9-point stencil
// https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1928&context=cstech
vec4 texLapl(sampler2D ch, vec2 uv) {
    return texStencil(ch, uv, float[](
        1.,   4., 1.,
        4., -20., 4.,
        1.,   4., 1.
    )) / 6.;
}

// horizontal gradient (Sobel filter)
vec4 texGradX(sampler2D ch, vec2 uv) {
    return texStencil(ch, uv, float[](
        -1., 0., 1.,
        -2., 0., 2.,
        -1., 0., 1.
    )) / 8.;
}

// vertical gradient (Sobel filter)
vec4 texGradY(sampler2D ch, vec2 uv) {
    return texStencil(ch, uv, float[](
         1.,  2.,  1.,
         0.,  0.,  0.,
        -1., -2., -1.
    )) / 8.;
}
```

## Buffer A

```glsl
vec4 char(vec2 p, int c) {
    if (p.x < 0. || p.x > 1. || p.y < 0.|| p.y > 1.) return vec4(0,0,0,1);
    return texture(iChannel2, p/16. + fract(vec2(c, 15-c/16)/16.));
}

void mainImage(out vec4 r, vec2 u) {
    vec2 i = u-u; r -= r;
    int z = 15;
    for(int k = (2*z+1)*(2*z+1); k-->0;) {
        i = vec2(k%(2*z+1),k/(2*z+1)) - float(z);
        float q = mix(.015, .06, smoothstep(0., 16e-3 * iResolution.y, iTime));
        r += .5 * texelFetch(iChannel0,ivec2(i+u),0) * (1. - q*dot(i,i))*exp(-q*dot(i,i));
    }
    vec2 uv = u/iResolution.xy;
    r = clamp(r,0.,1.);
    if (iFrame < 9) {
        if (u.x < .25 * iResolution.x) r.x = 1.;
        if (u.x > .75 * iResolution.x) r.y = 1.;
        //if (u.y < .05 * iResolution.y) r.z = 1.;
    }
    
    vec2 p = vec2(3.7,2.7) * (uv - vec2(.5,.38));
    float heart = distance(p, vec2(0, sqrt(abs(p.x))));
    if(heart > 1.) r = vec4(0);
}

```

## Image

```glsl
void mainImage( out vec4 r, in vec2 u )
{
    vec2 uv = u / iResolution.xy;
    
    // edge detection
    vec4 dx = texGradX(iChannel0, uv);
    vec4 dy = texGradY(iChannel0, uv);
    r = sqrt(dx*dx + dy*dy);
    
    // layering
    vec4 mask = texBlur(iChannel0, uv);
    float blend = smoothstep(.45, .55, uv.x);
    r.xy *= mask.xy - vec2(1. - blend, blend) * mask.yx;
    
    // colours
    r = 1. - clamp(2. * r, 0., 1.);
    r.y -= .5 - .5 * r.x;
    
    // paper
    r -= .05 * texture(iChannel1, .5 * u / iChannelResolution[1].xy).x;
    
    // shadow
    float shadow = .1 * length(texBlur(iChannel0, (u + vec2(-5, 5)) / iResolution.xy));
    shadow *= smoothstep(0., 4e-3 * iResolution.y, iTime);
    shadow *= 1. - clamp(mask.x + mask.y, 0., 1.);
    r -= shadow;
    
    //if (abs(u.x - .5 * iResolution.x) / iResolution.y > .5) r = vec4(0);
}
```
