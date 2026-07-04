
## Common

```glsl
//Orbit Noise
//by nimitz 2018 (twitter: @stormoid)

//See: https://www.shadertoy.com/view/Xt3cDn
uint baseHash(uvec2 p)
{
    uint h32 = p.y + 374761393U + p.x*3266489917U;
    h32 = 2246822519U*(h32^(h32 >> 15));
    h32 = 3266489917U*(h32^(h32 >> 13));
    return h32^(h32 >> 16);
}

uint baseHash(uvec3 p)
{
    p = 1103515245U*((p.xyz >> 1U)^(p.yzx));
    uint h32 = 1103515245U*((p.x^p.z)^(p.y>>3U));
    return h32^(h32 >> 16);
}

uint baseHash(uint p)
{
    p = 1103515245U*((p >> 1U)^(p));
    uint h32 = 1103515245U*((p)^(p>>3U));
    return h32^(h32 >> 16);
}

vec4 hash42(vec2 x)
{
    uint n = baseHash(floatBitsToUint(x));
    uvec4 rz = uvec4(n, n*16807U, n*n*48271U, n*n*69621U);
    return vec4(rz & uvec4(0x7fffffffU))/float(0x7fffffff);
}

vec2 hash21(float x)
{
    uint n = baseHash(uint(int(x)));
    uvec2 rz = uvec2(n, n*48271U); //see: http://random.mat.sbg.ac.at/results/karl/server/node4.html
    return vec2(rz.xy & uvec2(0x7fffffffU))/float(0x7fffffff);
}

```

## Image

```glsl
//Orbit Noise
//by nimitz 2018 (twitter: @stormoid)

/*
	This began while I was looking for ways to render high quality noise for
	offline noise pre-generation. I figured one of the main limitation with the
	standard noise algorithm was the rigid grid, which inevitably ends up creating
	linear artifacts in the produced noise.

	So a natural way of avoiding this issue is to let the underlying grid move
	somewhat freely. I started with experiments of moving cells on a grid and
	rendering them additively, with some extra code to improve color distribution.
	Then add an extra layer of variation by making the cells also have brightness values.
	This produces ok noise, but the underlying patterns were based on value noise,
	which is somewhat low quality.

	So to increase the quality, we can use two separate 2d vectors per cell,
	one for the cell displacement and the second to generate a (random) brightness gradient,
	this ends up creating noise that is pretty much gradient/coherent/perlin while the
	cell displacement is zeroed, which exhibits strong linear patterns. But with increased
	displacement, those artifacts disappear and the resulting noise is of high quality.
	As far as I'm aware, this technique is new so I'm taking the liberty of calling it "orbit noise"

	The noise works well with a 3 cell wide kernel (9 taps in 2d), but I am using
	a 5-wide kernel here just to help make the animation more satisfying. I am also including
	a 3D version of the noise to show that it is easy to extend to higher dimensions.
*/

//modified (simplified) nuttall filter for orbit noise
//see: https://www.shadertoy.com/view/XtVcWc
float nuttall(float x, float w)
{
    const float pi = 3.14159265358979;
    if (abs(x) > w)
        return 0.;
    //Standard Nuttall
    //return 0.355768 - 0.487396*cos(pi*x/w + pi) + 0.144232*cos(2.*pi*x/w) - 0.012604*cos(3.*pi*x/w + pi*3.);
    return 0.365 - (0.5)*cos(pi*x/w + pi) + (0.135)*cos(2.*pi*x/w);
}

float orbitNoise(vec2 p)
{
    vec2 ip = floor(p);
    vec2 fp = fract(p);
    float rz = 0.;
    float orbitRadius = .75;

    //16 taps
    for (int j = -1; j <= 2; j++)
    for (int i = -1; i <= 2; i++)
    {
        vec2 dp = vec2(j,i);
        vec4 rn = hash42(dp + ip) - 0.5;
        vec2 op = fp - dp + rn.zw*orbitRadius;
        rz += nuttall(length(op),1.85)*dot(rn.xy*1.7, op);
    }
    return rz*0.5+0.5;
}

//3D version, roughly 3 times as expensive as the 2d version
//still much cheaper than 3D voronoi since we can use a 3 cell wide kernel
float orbitNoise3D(vec3 p)
{
    vec3 ip = floor(p);
    vec3 fp = fract(p);
    float rz = 0.;
    float orbitRadius = 0.75; //Zero value for standard coherent/gradient/perlin noise

    for (int k = -1; k <= 2; k++)
    for (int j = -1; j <= 2; j++)
    for (int i = -1; i <= 2; i++)
    {
            vec3 dp = vec3(k,j,i);
            uint base = baseHash(floatBitsToUint(dp + ip));
        	vec3 rn1 = vec3(uvec3(base, base*16807U, base*48271U) & uvec3(0x7fffffffU))/float(0x7fffffff);
        	vec3 rn2 = vec3(base*1664525U, base*134775813U, base*22695477U) * (1.0/float(0xffffffffU)); //(2^32 LCGs)
        	vec3 op = fp - dp - (rn1.xyz-0.5)*orbitRadius;
        	rz += nuttall(length(op),1.85)*dot(rn2.xyz*1.0, op);
    }
    
    return rz*0.5 + 0.5;
}

//1D version, for completeness
float orbitNoise1D(float p)
{
    float ip = floor(p);
    float fp = fract(p);
    float rz = 0.;
    float orbitRadius = .75;

    for (int i = -1; i <= 2; i++)
    {
        vec2 rn = hash21(float(i) + ip) - 0.5;
        float op = fp - float(i) + rn.y*orbitRadius;
        rz += nuttall(abs(op), 1.85)*rn.x*3.0*op;
    }
    return rz*0.5+0.5;
}

float fbmNoise(vec2 p)
{
    float rz = 0.;
    float amp = 2.;
    for (int i = 0; i < 6; i++)
    {
		rz += orbitNoise(p)/amp;
        amp *= 2.;
        p *= 2.1;
        p += 12.5;
    }
    return rz;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 p = fragCoord/iResolution.xy;
    vec2 q = p;
	p.x *= iResolution.x/iResolution.y;
    
    p += iTime*.05;
    vec3 col = vec3(0);
    
    if (q.x < 0.53)
    	col = orbitNoise(p*36.3)*vec3(1);
        //col = orbitNoise3D(vec3(p*36.3,1.))*vec3(1);
    	//col = orbitNoise1D(p.x*50.)*vec3(1);
    else
    	col = fbmNoise(p*11.)*vec3(1);

    col *= smoothstep(0.,0.002, abs(q.x-0.53));

    fragColor = vec4(col,1.0);
}
```
