
## Common

```glsl
#define Bf(p) p
#define Bi(p) ivec2(p)
#define texel(a, p) texelFetch(a, Bi(p), 0)
#define pixel(a, p) texture(a, (p)/R)
#define ch0 iChannel0
#define ch1 iChannel1
#define ch2 iChannel2
#define ch3 iChannel3

#define PI 3.14159265

#define loop(i,x) for(int i = 0; i < x; i++)
#define range(i,a,b) for(int i = a; i <= b; i++)

#define dt 1.5

#define border_h 5.
vec2 R;
vec4 Mouse;
float time;

#define mass 1.

#define fluid_rho 0.5

float Pf(vec2 rho)
{
    //return 0.2*rho.x; //gas
    float GF = 1.;//smoothstep(0.49, 0.5, 1. - rho.y);
    return mix(0.5*rho.x,0.04*rho.x*(rho.x/fluid_rho - 1.), GF); //water pressure
}

mat2 Rot(float ang)
{
    return mat2(cos(ang), -sin(ang), sin(ang), cos(ang)); 
}

vec2 Dir(float ang)
{
    return vec2(cos(ang), sin(ang));
}


float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

float border(vec2 p)
{
    float bound = -sdBox(p - R*0.5, R*vec2(0.5, 0.5)); 
    float box = sdBox(Rot(0.*time)*(p - R*vec2(0.5, 0.6)) , R*vec2(0.05, 0.01));
    float drain = -sdBox(p - R*vec2(0.5, 0.7), R*vec2(1.5, 0.5));
    return max(drain,min(bound, box));
}

#define h 1.
vec3 bN(vec2 p)
{
    vec3 dx = vec3(-h,0,h);
    vec4 idx = vec4(-1./h, 0., 1./h, 0.25);
    vec3 r = idx.zyw*border(p + dx.zy)
           + idx.xyw*border(p + dx.xy)
           + idx.yzw*border(p + dx.yz)
           + idx.yxw*border(p + dx.yx);
    return vec3(normalize(r.xy), r.z + 1e-4);
}

uint pack(vec2 x)
{
    x = 65534.0*clamp(0.5*x+0.5, 0., 1.);
    return uint(round(x.x)) + 65535u*uint(round(x.y));
}

vec2 unpack(uint a)
{
    vec2 x = vec2(a%65535u, a/65535u);
    return clamp(x/65534.0, 0.,1.)*2.0 - 1.0;
}

vec2 decode(float x)
{
    uint X = floatBitsToUint(x);
    return unpack(X); 
}

float encode(vec2 x)
{
    uint X = pack(x);
    return uintBitsToFloat(X); 
}

struct particle
{
    vec2 X;
    vec2 V;
    vec2 M;
};
    
particle getParticle(vec4 data, vec2 pos)
{
    particle P; 
    P.X = decode(data.x) + pos;
    P.V = decode(data.y);
    P.M = data.zw;
    return P;
}

vec4 saveParticle(particle P, vec2 pos)
{
    P.X = clamp(P.X - pos, vec2(-0.5), vec2(0.5));
    return vec4(encode(P.X), encode(P.V), P.M);
}

vec3 hash32(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy+p3.yzz)*p3.zyx);
}

float G(vec2 x)
{
    return exp(-dot(x,x));
}

float G0(vec2 x)
{
    return exp(-length(x));
}

//diffusion amount
#define dif 1.12

vec3 distribution(vec2 x, vec2 p, float K)
{
    vec2 omin = clamp(x - K*0.5, p - 0.5, p + 0.5);
    vec2 omax = clamp(x + K*0.5, p - 0.5, p + 0.5); 
    return vec3(0.5*(omin + omax), (omax.x - omin.x)*(omax.y - omin.y)/(K*K));
}

/*
vec3 distribution(vec2 x, vec2 p, float K)
{
    vec4 aabb0 = vec4(p - 0.5, p + 0.5);
    vec4 aabb1 = vec4(x - K*0.5, x + K*0.5);
    vec4 aabbX = vec4(max(aabb0.xy, aabb1.xy), min(aabb0.zw, aabb1.zw));
    vec2 center = 0.5*(aabbX.xy + aabbX.zw); //center of mass
    vec2 size = max(aabbX.zw - aabbX.xy, 0.); //only positive
    float m = size.x*size.y/(K*K); //relative amount
    //if any of the dimensions are 0 then the mass is 0
    return vec3(center, m);
}*/

//diffusion and advection basically
void Reintegration(sampler2D ch, inout particle P, vec2 pos)
{
    //basically integral over all updated neighbor distributions
    //that fall inside of this pixel
    //this makes the tracking conservative
    range(i, -2, 2) range(j, -2, 2)
    {
        vec2 tpos = pos + vec2(i,j);
        vec4 data = texel(ch, tpos);
       
        particle P0 = getParticle(data, tpos);
       
        P0.X += P0.V*dt; //integrate position

        float difR = 0.9 + 0.21*smoothstep(fluid_rho*0., fluid_rho*0.333, P0.M.x);
        vec3 D = distribution(P0.X, pos, difR);
        //the deposited mass into this cell
        float m = P0.M.x*D.z;
        
        //add weighted by mass
        P.X += D.xy*m;
        P.V += P0.V*m;
        P.M.y += P0.M.y*m;
        
        //add mass
        P.M.x += m;
    }
    
    //normalization
    if(P.M.x != 0.)
    {
        P.X /= P.M.x;
        P.V /= P.M.x;
        P.M.y /= P.M.x;
    }
}

//force calculation and integration
void Simulation(sampler2D ch, inout particle P, vec2 pos)
{
    //Compute the SPH force
    vec2 F = vec2(0.);
    vec3 avgV = vec3(0.);
    range(i, -2, 2) range(j, -2, 2)
    {
        vec2 tpos = pos + vec2(i,j);
        vec4 data = texel(ch, tpos);
        particle P0 = getParticle(data, tpos);
        vec2 dx = P0.X - P.X;
        float avgP = 0.5*P0.M.x*(Pf(P.M) + Pf(P0.M)); 
        F -= 0.5*G(1.*dx)*avgP*dx;
        avgV += P0.M.x*G(1.*dx)*vec3(P0.V,1.);
    }
    avgV.xy /= avgV.z;

    //viscosity
    F += 0.*P.M.x*(avgV.xy - P.V);
    
    //gravity
    F += P.M.x*vec2(0., -0.0004);

     if(Mouse.z > 0.)
    {
        vec2 dm =(Mouse.xy - Mouse.zw*vec2(1.0,-1.0))/10.; 
        float d = distance(Mouse.xy, P.X)/20.;
        F += 0.001*dm*exp(-d*d);
       // P.M.y += 0.1*exp(-40.*d*d);
    }
    
    //integrate
    P.V += F*dt/P.M.x;

    //border 
    vec3 N = bN(P.X);
    float vdotN = step(N.z, border_h)*dot(-N.xy, P.V);
    P.V += 0.5*(N.xy*vdotN + N.xy*abs(vdotN));
    P.V += 0.*P.M.x*N.xy*step(abs(N.z), border_h)*exp(-N.z);
    
    if(N.z < 0.) P.V = vec2(0.);
    
    
    //velocity limit
    float v = length(P.V);
    P.V /= (v > 1.)?v:1.;
}



```

## Buffer A

```glsl

void mainImage( out vec4 U, in vec2 pos )
{
    R = iResolution.xy; time = iTime; Mouse = iMouse;
    ivec2 p = ivec2(pos);

    vec4 data = texel(ch0, pos); 
    
    particle P;// = getParticle(data, pos);
       
    Reintegration(ch0, P, pos);
   
    //initial condition
    if(iFrame < 1)
    {
        //random
        vec3 rand = hash32(pos);
        if(rand.z < 0.) 
        {
            P.X = pos;
            P.V = 0.5*(rand.xy-0.5) + vec2(0., 0.);
            P.M = vec2(mass, 0.);
        }
        else
        {
            P.X = pos;
            P.V = vec2(0.);
            P.M = vec2(1e-6);
        }
    }
    
    U = saveParticle(P, pos);
}
```


## Buffer B

```glsl
void mainImage( out vec4 U, in vec2 pos )
{
    R = iResolution.xy; time = iTime; Mouse = iMouse;
    ivec2 p = ivec2(pos);
        
    vec4 data = texel(ch0, pos); 
    
    particle P = getParticle(data, pos);
    
    
    if(P.M.x != 0.) //not vacuum
    {
        Simulation(ch0, P, pos);
    }
    
    if(length(P.X - R*vec2(0.8, 0.9)) < 10.) 
    {
        P.X = pos;
        P.V = 0.5*Dir(-PI*0.25 - PI*0.5 + 0.3*sin(0.4*time));
        P.M = mix(P.M, vec2(fluid_rho, 1.), 0.4);
    }

    if(length(P.X - R*vec2(0.2, 0.9)) < 10.) 
    {
        P.X = pos;
        P.V = 0.5*Dir(-PI*0.25 + 0.3*sin(0.3*time));
        P.M = mix(P.M, vec2(fluid_rho, 0.), 0.4);
    }
    
    U = saveParticle(P, pos);
}
```


## Buffer C

```glsl
//density

void mainImage( out vec4 fragColor, in vec2 pos )
{
    R = iResolution.xy; time = iTime;
    ivec2 p = ivec2(pos);

    vec4 data = texel(ch0, pos);
    particle P = getParticle(data, pos);
    
    //particle render
    vec4 rho = vec4(0.);
    range(i, -1, 1) range(j, -1, 1)
    {
        vec2 ij = vec2(i,j);
        vec4 data = texel(ch0, pos + ij);
        particle P0 = getParticle(data, pos + ij);

        vec2 x0 = P0.X; //update position
        //how much mass falls into this pixel
        rho += 1.*vec4(P.V, P.M)*G((pos - x0)/0.75); 
    }
    
    fragColor = rho;
}
```

## Image

```glsl
vec3 hsv2rgb( in vec3 c )
{
    vec3 rgb = clamp( abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),6.0)-3.0)-1.0, 0.0, 1.0 );

	rgb = rgb*rgb*(3.0-2.0*rgb); // cubic smoothing	

	return c.z * mix( vec3(1.0), rgb, c.y);
}

vec3 mixN(vec3 a, vec3 b, float k)
{
    return sqrt(mix(a*a, b*b, clamp(k,0.,1.)));
}

vec4 V(vec2 p)
{
    return pixel(ch1, p);
}

void mainImage( out vec4 col, in vec2 pos )
{
	R = iResolution.xy; time = iTime;
    //pos = R*0.5 + pos*0.1;
    ivec2 p = ivec2(pos);
    
    vec4 data = texel(ch0, pos);
    particle P = getParticle(data, pos);
    
    //border render
    vec3 Nb = bN(P.X);
    float bord = smoothstep(2.*border_h,border_h*0.5,border(pos));
    
    vec4 rho = V(pos);
    vec3 dx = vec3(-2., 0., 2.);
    vec4 grad = -0.5*vec4(V(pos + dx.zy).zw - V(pos + dx.xy).zw,
                         V(pos + dx.yz).zw - V(pos + dx.yx).zw);
    vec2 N = pow(length(grad.xz),0.2)*normalize(grad.xz+1e-5);
    float specular = pow(max(dot(N, Dir(1.4)), 0.), 3.5);
    float specularb = G(0.4*(Nb.zz - border_h))*pow(max(dot(Nb.xy, Dir(1.4)), 0.), 3.);
    
    float a = pow(smoothstep(fluid_rho*0., fluid_rho*2., rho.z),0.1);
    float b = exp(-1.7*smoothstep(fluid_rho*1., fluid_rho*7.5, rho.z));
    vec3 col0 = vec3(1., 0.5, 0.);
    vec3 col1 = vec3(0.1, 0.4, 1.);
	vec3 fcol = mixN(col0, col1, tanh(3.*(rho.w - 0.7))*0.5 + 0.5);
    // Output to screen
    col = vec4(3.);
    col.xyz = mixN(col.xyz, fcol.xyz*(1.5*b + specular*5.), a);
    col.xyz = mixN(col.xyz, 0.*vec3(0.5,0.5,1.), bord);
    col.xyz = tanh(col.xyz);
}
```
