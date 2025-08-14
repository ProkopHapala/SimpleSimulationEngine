#version 330 core
out vec4 fragColor;
uniform vec3 iResolution;
uniform int  iFrame;
uniform sampler2D iChannel0; // expects (vx,vy,density,vorticity)

vec3 palette(float t){
    // simple blue->cyan->yellow palette
    vec3 a=vec3(0.231,0.298,0.753);
    vec3 b=vec3(0.865,0.865,0.865);
    vec3 c=vec3(0.706,0.016,0.150);
    t=clamp(t,0.0,1.0);
    return mix(a, mix(b, c, smoothstep(0.4,1.0,t)), smoothstep(0.0,0.85,t));
}

void main(){
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    vec4 data = textureLod(iChannel0, uv, 0.0);
    float den = clamp(data.z/1.5, 0.0, 1.0);
    float vort = clamp(0.5 + 0.5*data.w*0.25, 0.0, 1.0);
    vec3 col = palette(den);
    col = mix(col, vec3(1.0,0.5,0.1), vort*0.15);
    //vec3 col = vec3(1.0,0.5,0.1);
    fragColor = vec4(col, 1.0);
}
