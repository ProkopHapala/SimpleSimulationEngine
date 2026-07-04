
## Common

```glsl

// https://nullprogram.com/blog/2018/07/31/
uint lowbias32 (uint x) {
    x = x ^ (x >> 16u);
    x = x * 0x7feb352du;
    x = x ^ (x >> 15u);
    x = x * 0x846ca68bu;
    x = x ^ (x >> 16u);
    return x;
}

uvec2 lowbias32(uvec2 v) {
    return uvec2(lowbias32(v.x), lowbias32(v.y));
}

// inverse
uint lowbias32_r (uint x) {
    x = x ^ (x >> 16u);
    x = x * 0x43021123u;
    x = x ^ (x >> 15u ^ x >> 30u);
    x = x * 0x1d69e2a5u;
    x = x ^ (x >> 16u);
    return x;
}

uvec2 lowbias32_r(uvec2 v) {
    return uvec2(lowbias32_r(v.x), lowbias32_r(v.y));
}

#define R2 19

#define SIGMA 1.414
#define M_PI 3.14159265359

float gaussian (float x, float sigma) {
    float h0 = x / sigma;
    float h = h0 * h0 * -0.5;
    float a = 1.0 / (sigma * sqrt(2.0 * M_PI));
    return a * exp(h);
}

float pow2(float x) {
    return x * x;
}

float distf(float v, float x) {
#if 0
    return 1.0 - abs(x);
#else
    return 1.0 / (1.0 + pow2(x));
#endif
}

struct line {
    float v[5];
};

vec2 quantify_error (sampler2D channel, ivec2 p, ivec2 sz, float val0, float val1) {
#if 1
    float Rf = float(R2) / 2.0;
    int R = int(Rf);
    float has0 = 0.0;
    float has1 = 0.0;
    float w = 0.0;

    //vec2 g = vec2(0.0, 1.0);
    for (int sy = -R; sy <= R; ++sy) {
        for (int sx = -R; sx <= R; ++sx) {
            float d = length(vec2(sx,sy));
            if ((d > Rf) || ((sx == 0) && (sy == 0)))
                continue;
            ivec2 t = (p + ivec2(sx,sy) + sz) % sz;            
			float v = texelFetch(channel, t, 0).r;
            //g += vec2(v, 1.0);

            float dist0 = (v - val0);
            float dist1 = (v - val1);

            float q = gaussian(d, SIGMA);

            w += q;            
            has0 += distf(val0, dist0) * q;
            has1 += distf(val1, dist1) * q;
            
        }
    }
    //vec2 avg = vec2(g.x + val0, g.x + val1) / g.y - 0.5;

    vec2 result = vec2(has0, has1) / w;
    //result = result * result;
    return result;
#else
    // FastNoise-inspired filter variance check
#define FETCH(OX, OY) texelFetch(channel, (p + ivec2((OX),(OY))) % sz, 0).r
    line lines[5];
    // cache
    for (int sy = -2; sy <= 2; ++sy) {
        for (int sx = -2; sx <= 2; ++sx) {
            lines[sy + 2].v[sx + 2] = FETCH(sx, sy);
        }
    }

    mat3 coeffs = mat3(
        1.0, 2.0, 1.0,
        2.0, 4.0, 2.0,
        1.0, 2.0, 1.0
    ) / 16.0;

    lines[2].v[2] = val0;
    float err[2] = float[](0.0, 0.0);
    for (int i = 0; i < 2; ++i) {
        // collect error of all neighboring filter kernels
        for (int sy = 0; sy < 3; ++sy) {
            for (int sx = 0; sx < 3; ++sx) {

                float acc = 0.0;
                // average
                for (int my = 0; my < 3; ++my) {
                    for (int mx = 0; mx < 3; ++mx) {
                        acc += lines[sy + my].v[sx + mx] * coeffs[my][mx];
                    }
                }           
                // variance
                float vr = 0.0;
                for (int my = 0; my < 3; ++my) {
                    for (int mx = 0; mx < 3; ++mx) {
                        float x = (lines[sy + my].v[sx + mx] - acc);
                        vr += x * x;
                    }
                }
                acc = (acc - 0.5) * (acc - 0.5);
                err[i] += acc + 1.0 / (1.0 + vr);
            }
        }
        lines[2].v[2] = val1;
    }
    return vec2(err[0], err[1]);
#endif
}

```

## Buffer A

```glsl

ivec2 flip(ivec2 p, uvec2 mask) {
    return ivec2(lowbias32_r(lowbias32(uvec2(p)) ^ mask));
}

float hash13(vec3 p3)
{
	p3  = fract(p3 * .1031);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    ivec2 sz = ivec2(iChannelResolution[0].xy);
    ivec2 p0 = ivec2(fragCoord);
    uvec2 mask = uvec2(lowbias32(uint(iFrame)));
    int M = 10 * 60;
    int F = (iFrame % M);
    float framef = float(F) / float(M);
    const float CHANCE_LIMIT = 0.618; // try to swap 62% of pixels
    if (F == 0) {
        int c = (p0.x * 61 + p0.y) % 256;
        fragColor = vec4(float(c) / 255.0, 0.0, 0.0, 1.0);
    } else {
        ivec2 p1 = flip(p0, mask);
        ivec2 pp0 = flip(p1, mask) % sz;
        p1 = p1 % sz;

        float chance0 = hash13(vec3(p0, float(iFrame)));
        float chance1 = hash13(vec3(p1, float(iFrame)));
        float chance = max(chance0, chance1);
        
        float v0 = texelFetch(iChannel0, p0, 0).r;
        float v1 = texelFetch(iChannel0, p1, 0).r;
        
        vec2 s0_x0 = quantify_error(iChannel0, p0, sz, v0, v1);
        vec2 s1_x1 = quantify_error(iChannel0, p1, sz, v1, v0);
        
        float err_s = s0_x0.x + s1_x1.x;
        float err_x = s0_x0.y + s1_x1.y;
        
        float p = v0;
        if ((chance < CHANCE_LIMIT) && (err_x < err_s)) {
            p = v1;
        }
        fragColor = vec4(p, 0.0, 0.0, 1.0);
    }
}                        

```

## Image

```glsl
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    ivec2 sz = ivec2(iResolution.xy);
    ivec2 px = ivec2(fragCoord);
   	float v_ref = texelFetch(iChannel1, ivec2(fragCoord)/1 % ivec2(iChannelResolution[1].xy), 0).r;
#if 0
    float v = texelFetch(iChannel0, px, 0).r;    
    float v_old = texelFetch(iChannel2, px, 0).r;
    if (uv.x > 0.5) {
        fragColor = vec4(v_ref,v_ref,v_ref,1.0);
    } else if (v != v_old) {
        fragColor = vec4(1.0,0.0,0.0,1.0);
    } else {
		fragColor = vec4(v,v,v,1.0);
    }
#else
    float v = texelFetch(iChannel0, px, 0).r;    
#if 0
    vec2 s0_x0 = quantify_error(iChannel0, px, sz, v, v);
    if (uv.x > 0.5) {
    v = s0_x0.x * 1.0;
    }
#endif

#if 0
    if (uv.x > 0.5) {
        v = v_ref; 
    }
    v = step(v, uv.y);
#endif
    //v = step(0.999/255.0, v);
    //v = step(v, 0.0);
    fragColor = vec4(v, v, v, 1.0);
#endif
}
```
