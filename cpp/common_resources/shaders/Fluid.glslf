

//Chimera's Breath
//by nimitz 2018 (twitter: @stormoid)
// https://www.shadertoy.com/view/4tGfDW
/*
	The main interest here is the addition of vorticity confinement with the curl stored in
	the alpha channel of the simulation texture (which was not used in the paper)
	this in turns allows for believable simulation of much lower viscosity fluids.
	Without vorticity confinement, the fluids that can be simulated are much more akin to
	thick oil.
	
	Base Simulation based on the 2011 paper: "Simple and fast fluids"
	(Martin Guay, Fabrice Colin, Richard Egli)
	(https://hal.inria.fr/inria-00596050/document)

	The actual simulation only requires one pass, Buffer A, B and C	are just copies 
	of each other to increase the simulation speed (3 simulation passes per frame)
	and Buffer D is drawing colors on the simulated fluid 
	(could be using particles instead in a real scenario)
*/

//#BEGIN_SHADER:SOLVER

#version 330 core


in       vec2      fUV;
uniform sampler2D  iChannel0; 
out vec4 gl_FragColor;

uniform float time;
uniform float dt;
uniform vec2  iResolution;
uniform float vorticity; 

// https://subscription.packtpub.com/book/game_development/9781782167020/1/ch01lvl1sec18/using-uniform-blocks-and-uniform-buffer-objects
// https://www.geeks3d.com/20140704/gpu-buffers-introduction-to-opengl-3-1-uniform-buffers-objects/
// https://paroj.github.io/gltut/Positioning/Tut07%20Shared%20Uniforms.html
/*
layout(std140) uniform ShaderToy{
  vec2  iResolution;
  float iTime;
  float iTimeStep;
  vec4      iFwColor;
  vec4      iBgColor;
  vec4      iMouse;
  sampler2D iChannel0;
  sampler2D iChannel1;
  sampler2D iChannel2;
  sampler2D iChannel3;
};
*/

#define USE_VORTICITY_CONFINEMENT
//#define MOUSE_ONLY

//Recommended values between 0.03 and 0.2
//higher values simulate lower viscosity fluids (think billowing smoke)
//#define VORTICITY_AMOUNT 0.11

float mag2(vec2 p){return dot(p,p);}
vec2 point1(float t) { t *= 0.12; return vec2(0.3,0.5 + sin(t         )*0.2); }
vec2 point2(float t) { t *= 0.22; return vec2(0.7,0.5 + cos(t + 1.5708)*0.2); }

vec4 solveFluid(sampler2D smp, vec2 uv, vec2 w ){
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
    
    
    vec2 newForce = vec2(0);
    newForce.xy += 2.75*vec2( .0003, 0.00005)/(mag2(uv-point1(time))+0.0001);
    newForce.xy += 2.75*vec2(-.0003, 0.00005)/(mag2(uv-point2(time))+0.0001);


    data.xy += dt*(viscForce.xy - K/dt*densDif + newForce); //update velocity
    //data.xy = max(vec2(0), abs(data.xy)-1e-4)*sign(data.xy); //linear velocity decay
    
    
    #ifdef USE_VORTICITY_CONFINEMENT
   	data.w    = (tr.y - tl.y - tu.x + td.x);
    vec2 vort = vec2(abs(tu.w) - abs(td.w), abs(tl.w) - abs(tr.w));
    vort     *= vorticity/length(vort + 1e-9)*data.w;
    data.xy  += vort;
    #endif
    
    
    data.y *= smoothstep(.5,.48,abs(uv.y-0.5)); //Boundaries
    
    data    = clamp(data, vec4(vec2(-10.), 0.5 , -10.), vec4(vec2(10.), 3.0 , 10.));

    //data = vec4(sin(uv*30.0),0.,1.);
    
    return data;
}

float length2(vec2 p){return dot(p,p);}
mat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,s,-s,c);}

void main( ){
    vec2 uv = fUV;
    vec2 w = 1.0/iResolution.xy;
    vec4 data = solveFluid(iChannel0, uv, w );
    gl_FragColor = data;
}








//#BEGIN_SHADER:RENDER

#version 330 core

in       vec2      fUV;
out vec4 gl_FragColor;

uniform float time;
uniform float dt;
uniform vec2  iResolution;
uniform sampler2D  iChannel0;
uniform sampler2D  iChannel1;

// https://paroj.github.io/gltut/Positioning/Tut07%20Shared%20Uniforms.html
/*
layout(std140) uniform ShaderToy{
  vec2  iResolution;
  float iTime;
  float idt;
  vec4      iFwColor;
  vec4      iBgColor;
  vec4      iMouse;
  sampler2D iChannel0;
  sampler2D iChannel1;
  sampler2D iChannel2;
  sampler2D iChannel3;
};
*/

mat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,s,-s,c);}

vec3 getPalette(float x, vec3 c1, vec3 c2, vec3 p1, vec3 p2){
    float x2 = fract(x/2.0);
    x = fract(x);   
    mat3 m = mat3(c1, p1, c2);
    mat3 m2 = mat3(c2, p2, c1);
    float omx = 1.0-x;
    vec3 pws = vec3(omx*omx, 2.0*omx*x, x*x);
    return clamp(mix(m*pws, m2*pws, step(x2,0.5)),0.,1.);
}

vec4 pal(float x){
    vec3 pal = getPalette(-x, vec3(0.2, 0.5, .7), vec3(.9, 0.4, 0.1), vec3(1., 1.2, .5), vec3(1., -0.4, -.0));
    return vec4(pal, 1.);
}

vec4 pal2(float x){
    vec3 pal = getPalette(-x, vec3(0.4, 0.3, .5), vec3(.9, 0.75, 0.4), vec3(.1, .8, 1.3), vec3(1.25, -0.1, .1));
    return vec4(pal, 1.);
}

//void mainImage( out vec4 fragColor, in vec2 fragCoord ){
void main( ){
    
    //vec2 uv = fragCoord.xy / iResolution.xy;
    //vec2 mo = iMouse.xy / iResolution.xy;
    //vec2 w = 1.0/iResolution.xy;
    
    vec2 uv = fUV;
    //vec2 w  = vec2(0.01,0.01);
    vec2 w = 1.0/iResolution.xy;
    
    vec2 velo = textureLod(iChannel0, uv            , 0.).xy;
    vec4 col  = textureLod(iChannel1, uv - dt*velo*w, 0.); //advection
    //vec4 col  = textureLod(iChannel1, uv - dt*vec2(1.0,0.5)*0.1, 0.); 
    
    //if (fragCoord.y < 1. && fragCoord.x < 1.)col = vec4(0);
    //vec4 lastMouse = texelFetch(iChannel1, ivec2(0,0), 0).xyzw;
    
    //if (iMouse.z > 1. && lastMouse.z > 1.){
    //    float str = smoothstep(-.5,1.,length(mo - lastMouse.xy/iResolution.xy));   
    //    col += str*0.0009/(pow(length(uv - mo),1.7)+0.002)*pal2(-iTime*0.7);
    //}
    //#ifndef MOUSE_ONLY
    //col += .0025/(0.0005+pow(length(uv - point1(iTime)),1.75))*dt*0.12*pal(iTime*0.05 - .0);
    //col += .0025/(0.0005+pow(length(uv - point2(iTime)),1.75))*dt*0.12*pal2(iTime*0.05 + 0.675);
    
    float bmix = 0.001;
    if(dot(col.xy,col.xy)<0.01) bmix=1.0;
    col.xy =  ( 1.0+sin(uv.yx*40.0) )*bmix + col.xy*(1.-bmix);
    //col.xy +=  (1.0+sin(uv.yx*40.0))*bmix;
    //#endif
    
    //if (iFrame < 20){  col = vec4(0.); }
    col = clamp(col, 0.,1.5);
    //col = max(col - (0.0001 + col*0.004)*.5, 0.); //decay
    //if (fragCoord.y < 1. && fragCoord.x < 1.) col = iMouse;
    gl_FragColor = col;
    //gl_FragColor = vec4( velo, 0.0, 1.0 );
    //gl_FragColor = vec4(sin(uv*30.),1.,1.);
    
    
    //gl_FragColor =  textureLod(iChannel0, fUV, 0.);
    //gl_FragColor =  textureLod(iChannel1, fUV, 0.);
}


