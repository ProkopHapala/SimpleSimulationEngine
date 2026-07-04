
## Common

```glsl
#define buf(p) texture(iChannel0,(p)/iResolution.xy)

#define N  vec2( 0, 1)
#define NE vec2( 1, 1)
#define E  vec2( 1, 0)
#define SE vec2( 1,-1)
#define S  vec2( 0,-1)
#define SW vec2(-1,-1)
#define W  vec2(-1, 0)
#define NW vec2(-1, 1)


/// from https://www.shadertoy.com/view/4tdSWr

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

const mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
float fbm(vec2 n) {
	float total = 0.0, amplitude = 0.1;
	for (int i = 0; i < 7; i++) {
		total += noise(n) * amplitude;
		n = m * n;
		amplitude *= 0.4;
	}
	return total;
}

///

```

## Buffer A

```glsl
// 2018 David A Roberts <https://davidar.io>

float slope(vec2 p, vec2 q) {
    return (buf(q).r - buf(p).r) / distance(p,q);
}

vec2 rec(vec2 p) { // direction of water flow at point
    vec2 d = N;
    if (slope(p + NE, p) > slope(p + d, p)) d = NE;
    if (slope(p + E,  p) > slope(p + d, p)) d = E;
    if (slope(p + SE, p) > slope(p + d, p)) d = SE;
    if (slope(p + S,  p) > slope(p + d, p)) d = S;
    if (slope(p + SW, p) > slope(p + d, p)) d = SW;
    if (slope(p + W,  p) > slope(p + d, p)) d = W;
    if (slope(p + NW, p) > slope(p + d, p)) d = NW;
    return d;
}

bool eq(vec2 p, vec2 q) {
    return distance(p,q) < 1e-3;
}

void mainImage( out vec4 r, in vec2 p ) {
    if (iFrame < 10 || iMouse.z > 0.) {
        r.r = clamp(5. * fbm(3. * p / iResolution.xy) + 0.5, 0., 1.);
        return;
    }
    r = buf(p);
    
    // flow accumulation
    r.g = 1.;
    if (eq(rec(p + N),  -N))  r.g += buf(p + N).g;
    if (eq(rec(p + NE), -NE)) r.g += buf(p + NE).g;
    if (eq(rec(p + E),  -E))  r.g += buf(p + E).g;
    if (eq(rec(p + SE), -SE)) r.g += buf(p + SE).g;
    if (eq(rec(p + S),  -S))  r.g += buf(p + S).g;
    if (eq(rec(p + SW), -SW)) r.g += buf(p + SW).g;
    if (eq(rec(p + W),  -W))  r.g += buf(p + W).g;
    if (eq(rec(p + NW), -NW)) r.g += buf(p + NW).g;
    
    // stream power
    vec4 receiver = buf(p + rec(p));
    float pslope = (r.r - receiver.r) / length(rec(p));
    r.r = max(r.r - 0.05 * pow(r.g, 0.8) * pow(pslope, 2.), receiver.r);
    
    // tectonic uplift
    r.r += 0.0004 * p.x/iResolution.x;
}

```

## Image

```glsl
#define PI 3.14159265359

void mainImage( out vec4 r, in vec2 p ) {
    float y = buf(p).r;
    vec2 grad = vec2(buf(p+E).r - buf(p+W).r, buf(p+N).r - buf(p+S).r);
    r = vec4(0.34, 0.52, 0.29, 1);
    r = mix(r, vec4(0.88, 0.85, 0.63, 1), smoothstep(0.500, 0.625, y));
	r = mix(r, vec4(0.93, 0.72, 0.40, 1), smoothstep(0.625, 0.750, y));
	r = mix(r, vec4(0.70, 0.60, 0.53, 1), smoothstep(0.750, 0.875, y));
	r = mix(r, vec4(1,1,1,1), smoothstep(0.875, 1., y));
    r.rgb *= 0.75 + 0.25 * cos(atan(grad.y, grad.x) + 0.25*PI);
    if (y < 0.5) r = mix(r, vec4(0.73, 0.80, 0.97, 1), 0.5);
}

```
