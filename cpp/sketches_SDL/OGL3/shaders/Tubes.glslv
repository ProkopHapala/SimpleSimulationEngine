attribute vec4  Position;
attribute vec3  Normal;
attribute float PathCoord;

varying   vec3  vPosition;
varying   vec3  vNormal;
varying  float  vPathCoord;

uniform   mat4  ModelviewProjection;

void main(){
    gl_Position = ModelviewProjection * Position;
    vPosition   = Position.xyz;
    vNormal     = Normal;
    vPathCoord  = PathCoord;
}
