
## Common

```glsl
// Noise simplex 2D by iq - https://www.shadertoy.com/view/Msf3WH

vec2 hash( vec2 p )
{
	p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p )
{
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

	vec2  i = floor( p + (p.x+p.y)*K1 );
    vec2  a = p - i + (i.x+i.y)*K2;
    float m = step(a.y,a.x); 
    vec2  o = vec2(m,1.0-m);
    vec2  b = a - o + K2;
	vec2  c = a - 1.0 + 2.0*K2;
    vec3  h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	vec3  n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
    return dot( n, vec3(70.0) );
}

```

## Buffer A

```glsl
// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

const int R = 10;         // space resolution = kernel radius
const float dt = .15;       // time step

const mat4 beta = mat4(.25,1,0,0, 1,.75,.75,0, 1,0,0,0, 0,0,0,0);   // kernel ring heights
const vec3 betaLen = vec3(2,3,1);    // kernel ring number
const vec3 mu = vec3(0.16, 0.22, 0.28);       // growth center
const vec3 sigma = vec3(0.025, 0.042, 0.025);   // growth width
const vec3 eta = vec3(1);     // growth strength

const float rho = 0.5;       // kernel center
const float omega = 0.15;    // kernel width

#define G(z) exp(-.5*(z)*(z))

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    ivec2 pos = ivec2(fragCoord);

    vec3 sum = vec3(0), total = vec3(0);
    for (int x = -R; x <= R; x++) for (int y = -R; y <= R; y++) {
        float r = length(vec2(x,y)/float(R));
        if (r > 1.) continue;
        float val = texelFetch(iChannel0, pos + ivec2(x,y), 0).x;
        vec3 height;
        for (int i = 0; i < 3; i++)
            height[i] = beta[i][int(r * betaLen[i])];
        vec3 weight = height * G((fract(r * betaLen) - rho) / omega);
        sum += val * weight;
        total += weight;
    }
    vec3 avg = sum / total;

    float r = texelFetch(iChannel0, pos, 0).x;
    r = mix(r, dot(G((avg - mu) / sigma), eta), dt);
    r = clamp(r, 0., 1.);

    if (iFrame < 2)
        r = noise(fragCoord/float(R) + iTime*100.);
    if (iMouse.z > 0.) {
        float d = length((fragCoord.xy - iMouse.xy) / iResolution.xx);
        if (d <= 5.*float(R)/iResolution.x)
        	r = noise(fragCoord/float(R) + iTime*100.);
    }

    fragColor = vec4(r);
}
```

## Image

```glsl
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = fragCoord.xy / iResolution.xy;
	fragColor = 1. - texture(iChannel0, uv);
}

```
