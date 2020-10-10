#version 330 core

in       vec2      fUV;
uniform  vec2      Const;    // julia constant
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

/*
    ShaderToy - https://www.shadertoy.com/view/XtcGD2
	See: A Simple Model of the Belousov-Zhabotinsky Reaction From First Principles
	by Alasdair Turner http://discovery.ucl.ac.uk/17241/1/17241.pdf
*/

#define TIMESTEP 0.01

//bool reset() {return texture(iChannel3, vec2(32.5/256.0, 0.5) ).x > 0.5;}

#define T(d) n += texture(texture_1, fract(vUv+d)).xyz;

//void mainImage( out vec4 fragColor, in vec2 fragCoord ) {

void main(){

    //vec2 vUv = fragCoord.xy / iResolution.xy;
    //vec4 t = vec4(1. / iResolution.xy, -1. / iResolution.y, 0.0);
    
    vec2 vUv = fUV;
    vec4 t = vec4(0.002, 0.002, -0.002, 0.0);
    
    vec3 p = texture(texture_1, vUv).xyz;
    vec3 n = vec3(0);
    
    // shorthand for summing the values over all 8 neighbors ... i.e. Diffusion
    T( t.wy) 
    T( t.xy) 
    T( t.xw) 
    T( t.xz) 
    T( t.wz) 
    T(-t.xy) 
    T(-t.xw) 
    T(-t.xz)
    
    // this line encodes the rules
    vec3 result = p + TIMESTEP * vec3(n.z - n.y, n.x - n.z, n.y - n.x);

    //if( p.xyz == vec3(0) ) {
    //    gl_FragColor = vec4( .5+sin(10000.0*fUV*fUV)*0.5, 0.0, 1.0);
    //}else{
        //gl_FragColor = vec4(clamp(result, 0.0, 1.0), 0.0);
        
        gl_FragColor = vec4(clamp((8.*p+n)/16.0, 0.0, 1.0), 0.0);
    //}
    
    //gl_FragColor = vec4( .5+sin(100.0*fUV)*0.5, 0.0, 1.0);

    // initialize with noise
    //if(p.xyz == vec3(0) || reset()) {
    //    fragColor = texture(iChannel1, vUv);
    //} else {
    //    fragColor = vec4(clamp(result, 0.0, 1.0), 0.0);
    //}
    
}

