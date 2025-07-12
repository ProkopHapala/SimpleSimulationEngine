//Chimera's Breath
//by nimitz 2018 (twitter: @stormoid)

//see "Common" tab for fluid simulation code

uniform vec3  iResolution;
// input textures
uniform sampler2D iChannel0;
uniform sampler2D iChannel1;

//#include "fluid/common.glslf"

mat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,s,-s,c);}

//shader incoming relating to this palette
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy / iResolution.xy;
    vec2 mo = iMouse.xy    / iResolution.xy;
    vec2 w  = 1.0/iResolution.xy;
    
    vec2 velo = textureLod(iChannel0, uv, 0.).xy;
    vec4 col  = textureLod(iChannel1, uv - dt*velo*w*3., 0.); //advection
    if (fragCoord.y < 1. && fragCoord.x < 1.) col = vec4(0);
    
    //vec4 lastMouse = texelFetch(iChannel1, ivec2(0,0), 0).xyzw;
    // if (iMouse.z > 1. && lastMouse.z > 1.)
    // {
    //     float str = smoothstep(-.5,1.,length(mo - lastMouse.xy/iResolution.xy));   
    //     col += str*0.0009/(pow(length(uv - mo),1.7)+0.002)*pal2(-iTime*0.7);
    // }
    // #ifndef MOUSE_ONLY
    // col += .0025/(0.0005+pow(length(uv - point1(iTime)),1.75))*dt*0.12*pal(iTime*0.05 - .0);
    // col += .0025/(0.0005+pow(length(uv - point2(iTime)),1.75))*dt*0.12*pal2(iTime*0.05 + 0.675);
    // #endif
    
    if (iFrame < 20){ col = vec4(0.); }
    
    col = clamp(col, 0.,5.);
    col = max(col - (0.0001 + col*0.004)*.5, 0.); //decay
    
    if (fragCoord.y < 1. && fragCoord.x < 1.) col = iMouse;

    fragColor = col;
    
}
