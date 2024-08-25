//
// Optimized signed distance for spherical inversion without fudge factor.
//
// Adjusts sd to be optimal, but can't be exact without assumptions about map().
//
// (c) timestamp @ shadertoy.com 2024-07
//

//
//          VIEWER DISCRETION IS ADVISED
//
//        Some spheres where turned inside 
//        out during this coding exercise.
//

//
//  Keyboard:
//
//    z    Toggle cubes vs spheres
//    x    Toggle 2D vs 3D
//    c    Show march step count in 3d
//
//  Mouse: sample point in 2D.
//

vec3 fold(vec3 p, float r, out vec2 v)
{
    v = vec2(dot(p, p), r * r);
    return v.y / v.x * p;    
}

float foldEnd(float sd, vec2 v)
{   
    return v.x / (v.y + sqrt(v.x) * abs(sd)) * sd;
}

//
//
// Demo
//
//


#define R                 iResolution.xy
#define Key_Z             0x5a
#define Key_X             0x58
#define Key_C             0x43

#define keyToggle(ascii)  (texelFetch(iChannel3,ivec2(ascii,2),0).x > 0.)
#define keyDown(ascii)    (texelFetch(iChannel3,ivec2(ascii,0),0).x > 0.)

bool doSpheres;

float map0(vec3 p)
{
    p.xy = mod(p.xy + .7, 1.4) - .7;
    vec3 b = abs(p) - .3;
    return doSpheres
        ? length(p) - 0.4
        : length(max(b,0.0)) + min(max(b.x,max(b.y,b.z)),0.0);    
}

float map1(vec3 p)
{
    vec2 v;
    vec3 q = fold(p, 1., v);
    float sd = map0(q);
    return foldEnd(sd, v);
}

float map(vec3 p)
{
    return max(map1(p * .5) / .5, abs(p.z) - 1.);
}

vec4 show2D(vec2 I)
{    
    vec2 p = (I+I-R)/R.y;    
    float h = 3.0 / iResolution.y;
    // derived from iq
    vec2 m = (2.0 * iMouse.xy - iResolution.xy) / iResolution.y;
    float d = map(vec3(p, 0));        
    vec3 col = 0.5 + 0.5 * sign(d) * (vec3(0.9, 0.6, 0.3) * 2. - 1.);    
    col *= 1.0 - exp(-20.0*abs(d));
    col *= 0.8 + 0.2*cos(140.0*d);
    col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,h,abs(d)) );
    if(iMouse.z < 0.001) m = vec2(cos(iTime * 0.4811), sin(iTime * 0.3211));
    d = map(vec3(m,0));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, h, abs(length(p-m)-abs(d))));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, h, length(p-m)-h));
    return vec4(col, 1);
}

vec3 getMapNormal(vec3 p)
{
    // iq
    vec2 e = vec2(1.0, -1.0) * 0.57735027 * .002;
    vec4 s = vec4(map(p + e.xyy)
        , map(p + e.yyx)
        , map(p + e.yxy)
        , map(p + e.xxx));
    return normalize((s.xzy + s.www) - (s.yxx + s.zyz));
}

vec2 rot(vec2 p, float a) { return vec2(p.x * cos(a) - p.y * sin(a), p.y * cos(a) + p.x * sin(a)); }

vec4 show3D(vec2 I)
{
    vec2 sc = (I+I-R)/R.y;    
    vec3 ro = vec3(-1.8,-0.6,1.7);
    ro.xy = rot(ro.xy, iTime * 0.4);
    vec3 rd = -normalize(ro);
    vec3 right = normalize(cross(rd, vec3(0, 0, 1)));
    vec3 up = normalize(cross(right, rd));
    rd = normalize(rd * 3. + right * sc.x + up * sc.y);

    int i, Iters = 100;
    float t = 0.;
    vec3 p, col;
    
    for(i=min(0, iFrame); i<Iters; i++)
    {
        p = ro + rd * t;
        float d = map(p);        
        if(d<.00001) break;
        t += d;
    }

    bool stp = keyToggle(Key_C);
    if(i < Iters)
    {    
        vec3 n = getMapNormal(p); 
        vec3 dif = vec3(0.945,0.820,0.475);
        vec3 lcol = vec3(1.);
        vec3 ld = normalize(vec3(1,0.5,2));    
        col = dif * lcol * max(0.005, dot(ld, n));
        if(stp) col.z = float(i) / float(Iters);
    }
    else if(stp && p.z > -2.)
        col.xz = vec2(1,0);
    
    return vec4(pow(col, vec3(1./2.2)), 1);
}

void mainImage( out vec4 O, in vec2 I)
{    
    doSpheres = keyToggle(Key_Z);
    O = keyToggle(Key_X) ? show3D(I) : show2D(I);    
}