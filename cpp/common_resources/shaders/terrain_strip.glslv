#version 330 core

layout(location = 0) in vec2 vUV;

noperspective out vec3 world_pos;
noperspective out vec3 world_nor;
smooth out vec3 fragColor;
//out float logz;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform vec3 camPos;
uniform mat4 camMat;

uniform vec2 uv0;
uniform vec2 pa;
uniform vec2 pb;
uniform vec3 mapScale;
uniform float derivScale;
uniform vec2  txStep;
uniform sampler2D txHeight;

vec4 bicubicSample( sampler2D tx, vec2 uv, vec2 d ){
    // see: https://stackoverflow.com/questions/20052381/glsl-performance-function-return-value-type
    vec4 p0 = texture2D( tx, uv );
    vec4 p1 = texture2D( tx, uv + vec2( d.x, d.y) );
    vec4 p2 = texture2D( tx, uv + vec2(-d.x, d.y) );
    vec4 p3 = texture2D( tx, uv + vec2( d.x,-d.y) );
    vec4 p4 = texture2D( tx, uv + vec2(-d.x,-d.y) );
    return (  2.0*p0  + p1 + p2 + p3 + p4)/6.0;
}

void main(){
    vec2 p      = pa*(1.0-vUV.x) + pb*vUV.x;
    //p        += normalize(p)*vUV.y + modelPos.xz;
    p           = p*vUV.y + modelPos.xz;
    //vec4 tx     = textureLod( txHeight, p*mapScale.xy + uv0, 0 );
    vec4 tx = bicubicSample( txHeight, p*mapScale.xy + uv0, txStep*0.125 );
    //world_nor   = tx.rgb;
    //world_nor   = (tx.rgb-0.5)*2.0;
    float deriv_sc = -2.0*mapScale.z*derivScale;
    world_nor   = normalize( vec3( (tx.x-0.5)*deriv_sc, 1.0, (tx.y-0.5)*deriv_sc  ) );
    
    float h     = tx.z * mapScale.z;
    world_pos   = vec3( p.x ,h+modelPos.y, p.y );
    
    //world_pos   = vec3( p.x ,0.0, p.y );
    //world_pos = vec3( vUV.x, 0.0, vUV.y );
    gl_Position = camMat * vec4( world_pos-camPos, 1 );
}
