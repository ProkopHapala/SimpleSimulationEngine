uniform vec3      iResolution;
uniform vec4      pos0;

vec3 gauss(vec2 uv, vec4 params){
    vec2  center = params.xy;
    float radius = params.z;
    float height = params.w;
    vec2  d = uv - center;
    float r2 = dot(d,d);
    float dGdx = -d.x/r2;
    float dGdy = -d.y/r2;
    float G = exp(-r2/(radius*radius));
    return vec3( dGdx,  1.,  dGdy )*height*G;
}

void main(){
    vec2 uv = gl_FragCoord.xy/iResolution.xy - 0.5;
    vec3 g = gauss(uv, pos0);
    g += gauss(uv, vec4(0.0,0.0,0.1,1.0));
    //g += +vec3(uv,0.0);
    gl_FragColor = vec4( g, 1. );
}
