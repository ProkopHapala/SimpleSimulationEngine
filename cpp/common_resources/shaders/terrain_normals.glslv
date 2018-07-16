#version 330 core

layout(location = 0) in vec3 vpos;

noperspective out vec3 world_pos;
noperspective out vec3 world_nor;
smooth out vec3 fragColor;
//out float logz;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform vec3 camPos;
uniform mat4 camMat;

uniform vec2  uv0;
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
    vec2  vUV = vpos.xy;
    float l_param = vpos.z;
    
    float muvx = (1.0-vUV.x);
    vec2 p     = (ps[0]*muvx + ps[1]*vUV.x)*(1.0-vUV.y) + (ps[2]*muvx + ps[3]*vUV.x)*vUV.y;
    p          += modelPos.xz;
    //vec4 tx   = textureLod( txHeight, p*mapScale.xy + uv0, 0 );
    //vec4 tx   = texture2D( txHeight, p*mapScale.xy + uv0, 0 );
    vec4 tx     = bicubicSample( txHeight, p*mapScale.xy + uv0, txStep*0.125 );

    vec2 deriv =  ( tx.xy-vec2(0.5) )*(-2.0*mapScale.z*derivScale);
    world_nor   = normalize( vec3( deriv.x, 1.0, deriv.y  ) );
    //p += deriv * flexibility * l_param;
    
    float h     = tx.z * mapScale.z;
    world_pos   = vec3( p.x, h+modelPos.y, p.y );
    
    world_pos += world_nor * flexibility * l_param;
    
    fragColor = vec3(1.0,0.0,0.0)*(1-l_param) +  vec3(0.0,0.0,1.0)*l_param;
    
    gl_Position = camMat * vec4( world_pos-camPos, 1 );
}
