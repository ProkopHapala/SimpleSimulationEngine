

## Buffer A

```glsl
/*
Author: Joe Eagar

https://www.shadertoy.com/view/wllSR4

pathtracer that has pixel decorrelation with blue noise
it's a very simple concept, see:

https://www.arnoldrenderer.com/research/dither_abstract.pdf

02/20/2020: added a preview mode.  tried to make navigation more interesting.

*/

//uncomment for preview mode
//#define PREVIEW

#ifdef PREVIEW
#define STEPS 50
#define IFSSTEPS 6
#define BOUNCES 2
#else
#define STEPS 100
#define IFSSTEPS 8
#define BOUNCES 4
#endif

#define COHERENCY_STEPS 19 //controls how may pixels in blue noise mask produce same seed
#define BLUE_NOISE_SEEDS 1
#define DSCALE 2.5 //scale per fractal level

#define SUN_ENERGY 1.0
#define SUN_SIZE 0.1
#define AMBIENT 1.0

//#define LIGHTS

struct State {
    vec3 origin;
    vec3 target;
    
    mat3 bases;
    float focal_dist;
};

State state;

struct Sample {
    float d;
    vec3 color;
    vec3 emission;
    vec3 no;
};

float tent(float f) {
    return 1.0 - abs(fract(f)-0.5)*2.0;
}

vec2 rot2d(vec2 v, float th) {
    return vec2(cos(th)*v.x + sin(th)*v.y, cos(th)*v.y - sin(th)*v.x);
}

vec3 cubenorm(vec3 co2) {
    vec3 aco2 = abs(co2);
            
    if (aco2[0] > aco2[1] && aco2[0] > aco2[2])
        co2 = vec3(sign(co2[0]), 0.0, 0.0);
    else if (aco2[1] > aco2[0] && aco2[1] > aco2[2])
        co2 = vec3(0.0, sign(co2[1]), 0.0);
        else
            co2 = vec3(0.0, 0.0, sign(co2[2]));
        
	return -co2;
}

float psuedoxor(float a, float b) {
	return a < 0.5 == b < 0.5 ? min(a, b) : max(a, b);
}

float cube(vec3 co2, float ll) {
    float l = max(abs(co2.x), abs(co2.y));
    l = max(l, abs(co2.z));
    
    return l - ll;
}

Sample s_cube(vec3 co, float l1) {
    Sample s;
    
    s.d = cube(co, l1);
    s.color = vec3(1.0,1.0,1.0)*0.78;
    s.no = -cubenorm(co);
    
    return s;
}


float sphere(vec3 co, float radius) {
    return length(co) - radius;
}

Sample s_sphere(vec3 co, float radius) {
    Sample s;
    
    s.d = sphere(co, radius);
    s.color = vec3(1.0, 1.0, 1.0)*0.78;
    s.no = normalize(co);
    
    return s;
}

void s_union(inout Sample a, Sample b) {
    if (a.d > b.d) {
        a = b;
    }
}

void s_diff(inout Sample a, Sample b) {
    if (a.d < -b.d) {
        a = b;
        a.d = -a.d;
        a.no = -a.no;
    }

}

void s_intersect(inout Sample a, Sample b) {
    if (a.d < b.d) {
        a = b;
    }

}


vec3 tent(vec3 f) {
    return vec3(tent(f.x), tent(f.y), tent(f.z));
}


float rand(float seed) {
    //seed += 50.0;
    
    float f = fract(sin(seed*3.14159265453)*59407.2751);
    
    f = fract(1.0 / (0.0000001 + 0.00001*f));
    
    return f;
}

vec3 randvec(float seed) {
    return vec3(
        rand(seed),
        rand(seed+12.23432),
        rand(seed+35.73423)
    ) - 0.5;
}

Sample fractal(vec3 co, float seed) {
    Sample s;
    
    vec3 center = vec3(1.0);
    vec3 co2 = co;
    int i;
    float th=0.0, thscale = 1.3, dscale = DSCALE;
    float k = 0.0, fi=0.0;

    //s = s_cube(co2, 0.5);
    s = s_sphere(co2, 0.5);
    
    float scale = dscale;
    
    for (i=0; i<IFSSTEPS; i++) {
        float x1, y1;
		
        /*
        co.xy = rot2d(co2.xy, th);
        co2.yz = rot2d(co2.yz, th);
        th += k;
        k += thscale*thscale;
		//*/
        
        co2 = floor(co*scale + 0.5)/scale;
        scale *= dscale;
        co2 = (co - co2);
        
        s_diff(s, s_cube(co2, 1.0/scale));
#ifdef LIGHTS
        float sz = 100.0;
        float f = tent(sz*co[0])*tent(sz*co[1])*tent(sz*co[2]);
        f = float(f > 0.75)*30.0;
        
        s.emission[1] = 0.25*f;
        s.emission[0] = f;
#endif       
        
    }
    
    return s;
}


Sample dsample(vec3 co, float seed) {
    return fractal(co, seed);
}


Sample trace(vec3 co, vec3 nray, out vec3 outco, out float found, float seed) {
    float dt = 3.0 / float(STEPS);
    Sample outl, f, l;
    
    found = 0.0;
    
    float facmul = pow(0.7, 1.0 / float(STEPS)), fac = 1.0;
    float t = 0.0;
    vec3 co2 = co;
    nray = normalize(nray);
    
    float mint=1000.0;
                     
    for (int i=0; i<STEPS; i++) {
        l = dsample(co2, seed);
        
        if (abs(l.d) < 0.0003) {// && abs(l.d) < mint) {//l < 0.1 && l > -0.1) { //abs(l) < 0.05) {
             outl = l;
		     outco = co2;
             found = 1.0;
             mint = abs(l.d);
             if (abs(l.d) < 0.00004)
            	break;
        }

        t += l.d; //*((cos(iTime*55.23532))*0.15+1.0);
        
        //don't let ray go behind origin
        t = max(t, 0.0);
        
        co2 = co + nray*t;
        fac *= facmul;
    }
    
    return outl;
}


/*
ambient occlusion.  currently unused.

samples the distance field at random points in a sphere
*/
float ao(vec3 co, vec3 nray, float max_dist, float seed) {
    int i;
    float tot=0.0;
    float sum=0.0;
    
    vec3 co2;
    
    for (int i=0; i<2; i++) {
        float seed2 = seed + 11.234 + 0.00123*float(i); //co.x*co.y*co.z + float(i)*3.14159;
        
        vec3 nray2 = normalize(randvec(seed2));
        float found;
        
        trace(co + nray2*0.005, nray2, co2, found, seed);
        sum += float(length(co2-co) < max_dist)*found;
        tot += 1.0;
    }
    
    sum = tot != 0.0 ? 1.0-sum / tot : 1.0;
    return sum;
}

//main pathtracing loop
vec3 shadeloop(vec3 co, vec3 nray, float seed) {
    vec3 co2;
    vec3 color;
    float seed2 = seed*5.0;
    
    for (int i=0; i<BOUNCES; i++) {
        seed2 += 1.32423;
        float found, f;
        Sample s = trace(co, nray, co2, found, seed2);
        float l = s.d;
		
        //sample ambient?
        if (found == 0.0) {
            color += vec3(0.4, 0.5, 0.9)*AMBIENT;
            break;
        }

        vec3 n = s.no;
		
        //planar area light instead of sun,
        //because it makes prettier shading :)
        vec3 light = vec3(3.0, 0.6, 3.0) + SUN_SIZE*randvec(seed2+1.0)*1.0*vec3(1.0,1.0,0.0); //nray[0], nray[1], nray[2]-2.5);
        vec3 ln = normalize(light-co2);
        
        f = max(dot(n, ln), 0.0)*SUN_ENERGY;
        
        vec3 co3;
        trace(co2 + -nray*0.004, ln, co3, found, seed2+2.0);

        f *= float(found < 0.1);
        
        color = color*s.color + s.color*f;
        color += s.emission;
        
        co = co2 - nray*0.004;
        nray = randvec(seed2+3.0);
        if (dot(nray, n) < 0.0) {
            nray = -nray;
        }
    }
    
    return color;
}

vec3 viewtrace(vec3 viewco, float seed) {
    float f = 0.0;

    vec3 planeco = state.origin + state.bases[2]*state.focal_dist;
    
    planeco += state.bases[0]*viewco[0];
    planeco += state.bases[1]*viewco[1];

    vec3 ray = planeco - state.origin;
    vec3 nray = normalize(ray);
    vec3 co2 = state.origin;// + nray*state.focal_dist + nray*1.2;
    vec3 co;
    
    return shadeloop(co2, nray, seed);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;
    vec4 inp = texture(iChannel0, uv);
    
    float seed = inp[3]*10.23423;
    
    //sample blue noise texture
    vec2 buv = fract(fragCoord.xy / iChannelResolution[1].xy);
    float bn = texture(iChannel1, buv)[0];
    bn = floor(bn*float(COHERENCY_STEPS))/float(COHERENCY_STEPS);
    
    float t = 1.0;
    vec2 viewco2d = fragCoord / iResolution.xy;
    viewco2d = (viewco2d - 0.5) * 2.0;
    viewco2d[0] *= iResolution.x / iResolution.y;
	
#if BLUE_NOISE_SEEDS
    //decorrelate pixel seeds with blue noise 
    seed += bn;
#else
    //white noise
    seed += rand(fragCoord.x*sqrt(3.0) + fragCoord.y*sqrt(5.0))*3.2432;
#endif 
    
    /*jitter the pixel grid a bit.
      for antiasing.  pseudo-uniform distribution (tent blur)
      by adding two random vectors.
     */
    viewco2d.xy += 2.0*(randvec(seed)+randvec(seed+1.1)).xy/iResolution.x;
    
    vec3 viewco = vec3(viewco2d, 0.0);
    vec3 clrout = vec3(0.0, 0.0, 0.0);
    float f = 0.0;
    
    viewco *= 1.5;
    
    //set camera origin
    float mx = iMouse.x/iResolution.x;
    float my = iMouse.y/iResolution.y;
    
    //offset so mouse at 0/0 is specific coordinates
    mx = fract(mx+0.73);
    //my = my+0.63;
    
    float th = -mx*1.5;
    float cr = 1.75 - 1.65*pow(abs(my), 0.2)*sign(my);
    
    float z = 2.0 - 1.0*(my-0.05);
    
    state.origin = vec3(cos(th), sin(th), z)*cr;
    state.target = vec3(0, 0, 0.33);
    
    state.focal_dist = 0.15;
    
    viewco[2] = state.focal_dist;
    viewco.xy *= 0.01;

    state.bases[2].xyz = normalize(state.target - state.origin);
    state.bases[0].xyz = normalize(cross(state.bases[2].xyz, vec3(0, 0, 1)));
    state.bases[1].xyz = normalize(cross(state.bases[0].xyz, state.bases[2].xyz));
    
	vec3 co = viewco;    
    float df = 1.0 / iResolution.y;

    clrout = viewtrace(co, seed);
    
    if (false || iMouse.z > 0.0) {
    	fragColor = vec4(clrout.rgb, 1.0);
    	return;
    }
    
    fragColor[3] = inp[3] + 1.0;
    fragColor.xyz = inp.xyz + clrout;
}
```

## Image

```glsl
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	vec2 uv = fragCoord / iResolution.xy;
    
    vec4 c;
    c = texture(iChannel0, uv);
    vec3 c2 = c.rgb/c.a;
    
    //toning
    c2 = sqrt(c2*2.0);
    //vec3 c3 = 1.0-exp(-c2*1.8);
    vec3 toning;
    //vec3 toning = mix(c2, c3, 0.5)*0.875;
    
    toning = exp(-(1.0-c2)*3.0);
    //toning = sqrt(c2)*0.5 + c2*c2*0.5;
    //toning = c.rgb/c.a;
    
    
    fragColor = vec4(toning, 1.0);
    
}
```