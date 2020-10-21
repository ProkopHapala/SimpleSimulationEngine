#version 330 core

in       vec2   fUV;

uniform  vec2       Const;    // julia constant
uniform  sampler2D  texture_1;

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


#define TIMESTEP 0.02        // TODO : This should be uniform parameter


out vec4 gl_FragColor;

/*
    ShaderToy - https://www.shadertoy.com/view/XtcGD2
	See: A Simple Model of the Belousov-Zhabotinsky Reaction From First Principles
	by Alasdair Turner http://discovery.ucl.ac.uk/17241/1/17241.pdf
*/




#define T(d) n += texture(texture_1, fract(vUv+d)).xyz;

void main(){

    vec2 vUv = fUV;
    vec4 t = vec4(0.001,0.0015, 0.0,0.0)*0.5;
    
    vec3 p = texture(texture_1, vUv).xyz;
    vec3 n = vec3(0);
    // shorthand for summing the values over all 8 neighbors ... i.e. Diffusion
    //T( t.wy)   T( t.wy)     T( t.xw)
    //T(-t.xw)                T( t.ww)
    //T(-t.xy)   T(-t.wy)     T(-t.xz)
               T( t.wy)     
    T(-t.xw)                T( t.xw)
               T(-t.wy)
    
    //  differential equation of B-Z reaction kintetics 
    vec3 result = p + TIMESTEP * vec3(n.z - n.y, n.x - n.z, n.y - n.x);

    gl_FragColor = vec4( result, 1.0 );

}

