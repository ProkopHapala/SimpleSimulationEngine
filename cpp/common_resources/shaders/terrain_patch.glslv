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

uniform vec2  uv0;
//uniform vec2  p00;
//uniform vec2  p01;
//uniform vec2  p10;
//uniform vec2  p11;
uniform vec2  ps[4];
uniform vec3  mapScale;
uniform float derivScale;
uniform vec2  txStep;

uniform sampler2D txHeight;

uniform float flexibility;


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
    float muvx = (1.0-vUV.x);
    //vec2 p     = (p00*muvx + p01*vUV.x)*(1.0-vUV.y) + (p10*muvx + p11*vUV.x)*vUV.y;
    vec2 p     = (ps[0]*muvx + ps[1]*vUV.x)*(1.0-vUV.y) + (ps[2]*muvx + ps[3]*vUV.x)*vUV.y;
    //p        += normalize(p)*vUV.y + modelPos.xz;
    p          += modelPos.xz;
    vec4 tx    = textureLod( txHeight, p*mapScale.xy + uv0, 0 );
    //vec4 tx = bicubicSample( txHeight, p*mapScale.xy + uv0, txStep*0.125 );

    //world_nor   = tx.rgb;
    //world_nor   = (tx.rgb-0.5)*2.0;
    
    vec2 deriv =  ( tx.xy-vec2(0.5) )*(-2.0*mapScale.z*derivScale);
    world_nor   = normalize( vec3( deriv.x, 1.0, deriv.y  ) );
    //world_nor   = vec3( tx.x, 0.0, tx.y );
    //world_nor   = vec3( 0.0,1.0,0.0 );
    
    p += deriv * flexibility;   // makes hil-tops thinner and sharper
    
    
    float h     = tx.z * mapScale.z;
    world_pos   = vec3( p.x, h+modelPos.y, p.y );
    
    //world_pos   = vec3( p.x ,0.0, p.y );
    //world_pos = vec3( vUV.x, 0.0, vUV.y );
    gl_Position = camMat * vec4( world_pos-camPos, 1 );
}
