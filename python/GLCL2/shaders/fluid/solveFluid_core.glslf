#version 330 core
out vec4 fragColor;
uniform vec3  iResolution;
uniform int   iFrame;
uniform sampler2D iChannel0;
uniform vec4 driver;

//#include "fluid/common.glslf"

float length2(vec2 p){return dot(p,p);}
mat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,s,-s,c);}

//Chimera's Breath
//by nimitz 2018 (twitter: @stormoid)
/*
    The main interest here is the addition of vorticity confinement with the curl stored in
    the alpha channel of the simulation texture (which was not used in the paper)
    this in turns allows for believable simulation of much lower viscosity fluids.
    Without vorticity confinement, the fluids that can be simulated are much more akin to
    thick oil.
    
    Base Simulation based on the 2011 paper: "Simple and fast fluids"
    (Martin Guay, Fabrice Colin, Richard Egli)
    (https://hal.inria.fr/inria-00596050/document)

    The actual simulation only requires one pass, Buffer A, B and C are just copies 
    of each other to increase the simulation speed (3 simulation passes per frame)
    and Buffer D is drawing colors on the simulated fluid 
    (could be using particles instead in a real scenario)
*/

#define dt 0.15
#define USE_VORTICITY_CONFINEMENT
//#define MOUSE_ONLY

//Recommended values between 0.03 and 0.2
//higher values simulate lower viscosity fluids (think billowing smoke)
#define VORTICITY_AMOUNT 0.11

float mag2(vec2 p){return dot(p,p);} 

vec4 solveFluid(sampler2D smp, vec2 uv, vec2 w ){
    // texture components are (vx,vy,density,vorticity)
    const float K = 0.2;
    const float v = 0.55;

    vec4 data = textureLod(smp, uv, 0.0);
    vec4 tr = textureLod(smp, uv + vec2(w.x , 0), 0.0);
    vec4 tl = textureLod(smp, uv - vec2(w.x , 0), 0.0);
    vec4 tu = textureLod(smp, uv + vec2(0 , w.y), 0.0);
    vec4 td = textureLod(smp, uv - vec2(0 , w.y), 0.0);
    
    vec3 dx = (tr.xyz - tl.xyz)*0.5;
    vec3 dy = (tu.xyz - td.xyz)*0.5;
    vec2 densDif = vec2(dx.z ,dy.z);
    
    data.z -= dt*dot(vec3(densDif, dx.x + dy.y) ,data.xyz); //density
    vec2 laplacian = tu.xy + td.xy + tr.xy + tl.xy - 4.0*data.xy;
    vec2 viscForce = vec2(v)*laplacian;
    data.xyw = textureLod(smp, uv - dt*data.xy*w, 0.).xyw; //advection

    vec2 newForce = vec2(0.);
    newForce += 0.00001*driver.zw/(mag2(uv-driver.xy)+0.0001);

    data.xy += dt*(viscForce.xy - K/dt*densDif + newForce); //update velocity
    data.xy = max(vec2(0), abs(data.xy)-1e-4)*sign(data.xy); //linear velocity decay
#ifdef USE_VORTICITY_CONFINEMENT
    data.w = (tr.y - tl.y - tu.x + td.x);
    vec2 vort = vec2(abs(tu.w) - abs(td.w), abs(tl.w) - abs(tr.w));
    vort *= VORTICITY_AMOUNT/length(vort + 1e-9)*data.w;
    data.xy += vort;
#endif
    data.y *= smoothstep(.5,.48,abs(uv.y-0.5)); //Boundaries
    data = clamp(data, vec4(vec2(-10), 0.5 , -10.), vec4(vec2(10), 3.0 , 10.));

    return data;
}

void main(){
    vec2 uv = gl_FragCoord.xy/iResolution.xy;
    vec2 w  = 1.0/iResolution.xy;
    vec4 data;
    if(iFrame<20){ data = vec4(0.5,0.0,1.0,0.0); }
    else         { data = solveFluid(iChannel0, uv, w ); }
    fragColor = data;
}
