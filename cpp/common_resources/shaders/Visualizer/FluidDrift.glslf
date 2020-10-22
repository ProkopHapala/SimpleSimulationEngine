

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

//#BEGIN_SHADER:RENDER

#version 330 core

in       vec2      fUV;
out vec4 gl_FragColor;

uniform float iTime;
uniform float iTimeStep;
uniform vec2  iResolution;
uniform sampler2D  iChannel0; 
uniform sampler2D  iChannel1; 
uniform sampler2D  iChannel2; 
uniform sampler2D  iChannel3; 

// https://paroj.github.io/gltut/Positioning/Tut07%20Shared%20Uniforms.html
/*
layout(std140) uniform ShaderToy{
  vec2  iResolution;
  float iTime;
  float iiTimeStep;
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
    
    vec2 velo = textureLod(iChannel1, uv                   , 0.).xy;
    vec4 col  = textureLod(iChannel0, uv - iTimeStep*velo*w, 0.); //advection
    //vec4 col  = textureLod(iChannel0, uv - iTimeStep*vec2(1.0,0.5)*0.1, 0.); 
    
    //if (fragCoord.y < 1. && fragCoord.x < 1.)col = vec4(0);
    //vec4 lastMouse = texelFetch(iChannel1, ivec2(0,0), 0).xyzw;
    
    //if (iMouse.z > 1. && lastMouse.z > 1.){
    //    float str = smoothstep(-.5,1.,length(mo - lastMouse.xy/iResolution.xy));   
    //    col += str*0.0009/(pow(length(uv - mo),1.7)+0.002)*pal2(-iTime*0.7);
    //}
    //#ifndef MOUSE_ONLY
    //col += .0025/(0.0005+pow(length(uv - point1(iTime)),1.75))*iTimeStep*0.12*pal(iTime*0.05 - .0);
    //col += .0025/(0.0005+pow(length(uv - point2(iTime)),1.75))*iTimeStep*0.12*pal2(iTime*0.05 + 0.675);
    
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


