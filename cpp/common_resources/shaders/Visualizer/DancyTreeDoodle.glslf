

// https://www.shadertoy.com/view/wslGz7

#version 330 core


in       vec2      fUV;
uniform sampler2D  iChannel0; 
out vec4 gl_FragColor;

uniform float iTime;
//uniform float dt;
uniform vec2  iResolution;

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

vec2 po (vec2 v) {
	return vec2(length(v),atan(v.y,v.x));
}

vec2 ca (vec2 u) {
	return u.x*vec2(cos(u.y),sin(u.y));
}

float ln (vec2 p, vec2 a, vec2 b) { 
    float r = dot(p-a,b-a)/dot(b-a,b-a);
    r = clamp(r,0.,1.);
    p.x+=(0.7+0.5*sin(0.1*iTime))*0.2*smoothstep(1.,0.,abs(r*2.-1.))*sin(3.14159*(r-4.*iTime));
    return (1.+0.5*r)*length(p-a-(b-a)*r);
}

//void mainImage( out vec4 Q, in vec2 U )
void main(  ){   
    vec4 Q; 
    vec2 U = fUV*4.0-2.0;
    vec2 R = iResolution.xy;
    U.x*=(iResolution.x/iResolution.y);
 	float r = 1e9;
 	//U = 4.*(U-0.5*R)/R.y;
 	U.y += 1.5;
 	Q = vec4(0);
 	for (int i = 1; i < 20; i++) {
        U = ca( po(U) + 0.3*( sin(2.*iTime) + 0.5*sin(4.53*iTime) + 0.1*cos(12.2*iTime))*vec2(0,1) );
        r = min(r,ln(U,vec2(0),vec2(0,1.)));
        U.y-=1.;
        
        U.x  =  abs(U.x);
        U   *=  1.4+0.1*sin(iTime) + 0.05*sin(0.2455*iTime)*(float(i));
        U    =  po(U);
        U.y +=  1. + 0.5*sin(0.553*iTime)*sin(sin(iTime)*float(i) ) + 0.1*sin(0.4*iTime) + 0.05*sin(0.554*iTime);
        U    =  ca(U);
        
        Q+=sin(1.5*exp(-1e2*r*r)*1.4*vec4(1,-1.8,1.9,4)+iTime);
        
 	}
 	Q/=18.;

    gl_FragColor = Q;
}

