
## Common

```glsl
#define buf(p) texture(iChannel0,(p)/iResolution.xy)

#define N vec2( 0, 1)
#define E vec2( 1, 0)
#define S vec2( 0,-1)
#define W vec2(-1, 0)

#define PI 3.14159265359


// https://www.shadertoy.com/view/4tdSWr
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License

const mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );

vec2 hash( vec2 p ) {
	p = vec2(dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)));
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p ) {
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;
	vec2 i = floor(p + (p.x+p.y)*K1);	
    vec2 a = p - i + (i.x+i.y)*K2;
    vec2 o = (a.x>a.y) ? vec2(1.0,0.0) : vec2(0.0,1.0); //vec2 of = 0.5 + 0.5*vec2(sign(a.x-a.y), sign(a.y-a.x));
    vec2 b = a - o + K2;
	vec2 c = a - 1.0 + 2.0*K2;
    vec3 h = max(0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
    return dot(n, vec3(70.0));	
}

float fbm(vec2 n) {
	float total = 0.0, amplitude = 0.1;
	for (int i = 0; i < 7; i++) {
		total += noise(n) * amplitude;
		n = m * n;
		amplitude *= 0.4;
	}
	return total;
}


// Hash without Sine
// Creative Commons Attribution-ShareAlike 4.0 International Public License
// Created by David Hoskins.

// https://www.shadertoy.com/view/4djSRW
// Trying to find a Hash function that is the same on ALL systens
// and doesn't rely on trigonometry functions that change accuracy 
// depending on GPU. 
// New one on the left, sine function on the right.
// It appears to be the same speed, but I suppose that depends.

// * Note. It still goes wrong eventually!
// * Try full-screen paused to see details.


#define ITERATIONS 4


// *** Change these to suit your range of random numbers..

// *** Use this for integer stepped ranges, ie Value-Noise/Perlin noise functions.
#define HASHSCALE1 .1031
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define HASHSCALE4 vec4(.1031, .1030, .0973, .1099)

// For smaller input rangers like audio tick or 0-1 UVs use these...
//#define HASHSCALE1 443.8975
//#define HASHSCALE3 vec3(443.897, 441.423, 437.195)
//#define HASHSCALE4 vec3(443.897, 441.423, 437.195, 444.129)



//----------------------------------------------------------------------------------------
//  1 out, 1 in...
float hash11(float p)
{
	vec3 p3  = fract(vec3(p) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

//----------------------------------------------------------------------------------------
//  1 out, 2 in...
float hash12(vec2 p)
{
	vec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

//----------------------------------------------------------------------------------------
//  1 out, 3 in...
float hash13(vec3 p3)
{
	p3  = fract(p3 * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

//----------------------------------------------------------------------------------------
//  2 out, 1 in...
vec2 hash21(float p)
{
	vec3 p3 = fract(vec3(p) * HASHSCALE3);
	p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.xx+p3.yz)*p3.zy);

}

//----------------------------------------------------------------------------------------
///  2 out, 2 in...
vec2 hash22(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);

}

//----------------------------------------------------------------------------------------
///  2 out, 3 in...
vec2 hash23(vec3 p3)
{
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}

//----------------------------------------------------------------------------------------
//  3 out, 1 in...
vec3 hash31(float p)
{
   vec3 p3 = fract(vec3(p) * HASHSCALE3);
   p3 += dot(p3, p3.yzx+19.19);
   return fract((p3.xxy+p3.yzz)*p3.zyx); 
}


//----------------------------------------------------------------------------------------
///  3 out, 2 in...
vec3 hash32(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}

//----------------------------------------------------------------------------------------
///  3 out, 3 in...
vec3 hash33(vec3 p3)
{
	p3 = fract(p3 * HASHSCALE3);
    p3 += dot(p3, p3.yxz+19.19);
    return fract((p3.xxy + p3.yxx)*p3.zyx);

}

//----------------------------------------------------------------------------------------
// 4 out, 1 in...
vec4 hash41(float p)
{
	vec4 p4 = fract(vec4(p) * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
    
}

//----------------------------------------------------------------------------------------
// 4 out, 2 in...
vec4 hash42(vec2 p)
{
	vec4 p4 = fract(vec4(p.xyxy) * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);

}

//----------------------------------------------------------------------------------------
// 4 out, 3 in...
vec4 hash43(vec3 p)
{
	vec4 p4 = fract(vec4(p.xyzx)  * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

//----------------------------------------------------------------------------------------
// 4 out, 4 in...
vec4 hash44(vec4 p4)
{
	p4 = fract(p4  * HASHSCALE4);
    p4 += dot(p4, p4.wzxy+19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

```

## Buffer A

```glsl
// 2018 David A Roberts <https://davidar.io>

// plate movement
vec2 move(vec2 v, int iFrame) {
    if (hash11(float(iFrame)) < 0.2 && hash13(vec3(v,iFrame)) < length(v)) {
        if (hash13(vec3(v,iFrame)) < abs(v.x) / (abs(v.x) + abs(v.y))) {
            return vec2(sign(v.x),0.);
        } else {
            return vec2(0.,sign(v.y));
        }
    }
    return vec2(0);
}
#define MOVE(c) move(c.xy, iFrame)

void mainImage(out vec4 c, in vec2 p) {
    if(iFrame < 10 || iMouse.z > 0.) {
        c = vec4(0);
        c.z = 15. * clamp(5. * fbm(3. * p / iResolution.xy) + 0.5, 0., 1.);
        return;
    }
    
    c = buf(p);
    vec4 n = buf(p + N);
    vec4 e = buf(p + E);
    vec4 s = buf(p + S);
    vec4 w = buf(p + W);
    
    // erosion
    c.z = max(0., c.z + 0.01 * (e.z + w.z + n.z + s.z - 4.*c.z));
    
    c.z += 0.0015;
    
    if(iFrame % 500 < 10) {
        // generate new plate boundaries
        c.xy = vec2(0);
    } else if(c.xy == vec2(0)) { // no plate under this point yet
        if(length(hash33(vec3(p,iFrame))) < 0.01) {
            // seed a new plate with random velocity
            c.xy = hash23(vec3(p,iFrame)) - 0.5;
        } else {
            // accretion
            int dir = int(4.*hash13(vec3(p,iFrame)));
            if(dir == 0) c.xy = s.xy;
            if(dir == 1) c.xy = w.xy;
            if(dir == 2) c.xy = n.xy;
            if(dir == 3) c.xy = e.xy;
        }
    } else if (MOVE(n) == S) {
        if (MOVE(c) != S) n.z += 1. + c.z; // subduction
        c = n;
    } else if (MOVE(e) == W) {
        if (MOVE(c) != W) e.z += 1. + c.z; // subduction
        c = e;
    } else if (MOVE(s) == N) {
        if (MOVE(c) != N) s.z += 1. + c.z; // subduction
        c = s;
    } else if (MOVE(w) == E) {
        if (MOVE(c) != E) w.z += 1. + c.z; // subduction
        c = w;
    } else if (MOVE(c) != vec2(0) && buf(p - MOVE(c)).xy != vec2(0)) { 
        // rift
        c = vec4(0);
    }
}
```

## Image

```glsl
#define DEBUG false

void mainImage(out vec4 r, in vec2 p) {
    float y = buf(p).z / 50.;
    if (DEBUG) {
        vec4 c = buf(p);
        r.rgb = (c.xy == vec2(0)) ? vec3(1) : .6 + .6 * cos(atan(c.y,c.x) + vec3(0,23,21));
        r.rgb *= 0.1 + 0.9 * y;
    } else if (y < 0.15) { // ocean
        r = mix(vec4(0.01, 0.02, 0.08, 1), vec4(0.11, 0.28, 0.51, 1), y/0.15);
    } else { // land
        vec2 grad = vec2(buf(p+E).z - buf(p+W).z, buf(p+N).z - buf(p+S).z);
        r = vec4(0.08, 0.14, 0.03, 1);
        r = mix(r, vec4(0.18, 0.26, 0.08, 1), smoothstep(0.15, 0.25, y));
        r = mix(r, vec4(0.52, 0.39, 0.26, 1), smoothstep(0.25, 0.50, y));
        r = mix(r, vec4(0.32, 0.3, 0.2, 1), smoothstep(0.50, 0.75, y));
        r = mix(r, vec4(1,1,1,1), smoothstep(0.75, 1., y));
        r.rgb *= 0.75 + 0.25 * cos(atan(grad.y, grad.x) + 0.25*PI)
                             * clamp(0.2 * length(grad), 0., 1.);
    }
}
```
