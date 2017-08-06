#extension GL_EXT_geometry_shader4 : enable

varying in  vec3 vPosition[4];
varying in  vec3 vNormal  [4];

varying out vec3 gPosition;
varying out vec3 gEndpoints[4];
varying out vec3 gEndplanes[2];

uniform float Radius;
uniform mat4  Modelview;
uniform mat4  Projection;

vec4 obb     [8];
vec4 obbPrime[8];

bool isFront(int a, int b, int c){
    vec3 i = vec3(obbPrime[b].xy - obbPrime[a].xy, 0);
    vec3 j = vec3(obbPrime[c].xy - obbPrime[a].xy, 0);
    return cross(i, j).z > 0.0;
}

void emit(int a, int b, int c, int d){
    gPosition = obb[a].xyz; gl_Position = obbPrime[a]; EmitVertex();
    gPosition = obb[b].xyz; gl_Position = obbPrime[b]; EmitVertex();
    gPosition = obb[c].xyz; gl_Position = obbPrime[c]; EmitVertex();
    gPosition = obb[d].xyz; gl_Position = obbPrime[d]; EmitVertex();
}

void main(){
    // Pass raytracing inputs to fragment shader:
    vec3 p0, p1, p2, p3, n0, n1, n2;
    p0 = (Modelview * vec4(vPosition[0], 1)).xyz;
    p1 = (Modelview * vec4(vPosition[1], 1)).xyz;
    p2 = (Modelview * vec4(vPosition[2], 1)).xyz;
    p3 = (Modelview * vec4(vPosition[3], 1)).xyz;
    n0 = normalize(p1-p0);
    n1 = normalize(p2-p1);
    n2 = normalize(p3-p2);
    gEndpoints[0] = p0; gEndpoints[1] = p1;
    gEndpoints[2] = p2; gEndpoints[3] = p3;
    gEndplanes[0] = normalize(n0+n1);
    gEndplanes[1] = normalize(n1+n2);

    // Compute object-space plane normals:
    p0 = vPosition[0]; p1 = vPosition[1];
    p2 = vPosition[2]; p3 = vPosition[3];
    n0 = normalize(p1-p0);
    n1 = normalize(p2-p1);
    n2 = normalize(p3-p2);
    vec3 u = normalize(n0+n1);
    vec3 v = normalize(n1+n2);

    // Generate a basis for the cuboid:
    vec3 j = n1;
    vec3 i = vNormal[1];
    vec3 k = cross(i, j);

    // Compute the eight corners:
    float r = Radius; float d;
    d = 1.0/dot(u,j); p1 -= j*r*sqrt(d*d-1.0);
    d = 1.0/dot(v,j); p2 += j*r*sqrt(d*d-1.0);
    obb[0] = Modelview*vec4(p1 + i*r + k*r,1);
    obb[1] = Modelview*vec4(p1 + i*r - k*r,1);
    obb[2] = Modelview*vec4(p1 - i*r - k*r,1);
    obb[3] = Modelview*vec4(p1 - i*r + k*r,1);
    obb[4] = Modelview*vec4(p2 + i*r + k*r,1);
    obb[5] = Modelview*vec4(p2 + i*r - k*r,1);
    obb[6] = Modelview*vec4(p2 - i*r - k*r,1);
    obb[7] = Modelview*vec4(p2 - i*r + k*r,1);
    for (int i = 0; i < 8; i++) obbPrime[i] = Projection * obb[i];
    
    // Emit the front faces of the cuboid:
    if (isFront(0,1,3)) emit(0,1,3,2); EndPrimitive();
    if (isFront(5,4,6)) emit(5,4,6,7); EndPrimitive();
    if (isFront(4,5,0)) emit(4,5,0,1); EndPrimitive();
    if (isFront(3,2,7)) emit(3,2,7,6); EndPrimitive();
    if (isFront(0,3,4)) emit(0,3,4,7); EndPrimitive();
    if (isFront(2,1,6)) emit(2,1,6,5); EndPrimitive();
}

