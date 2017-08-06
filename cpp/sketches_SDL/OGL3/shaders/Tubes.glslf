uniform vec4 Color;

varying vec3 gEndpoints[4];
varying vec3 gEndplanes[2];
varying vec3 gPosition;

uniform float Radius;
uniform mat4  Projection;
uniform vec3  LightDirection;
uniform vec3  DiffuseMaterial;
uniform vec3  AmbientMaterial;
uniform vec3  SpecularMaterial;
uniform float Shininess;

vec3 perp(vec3 v){
    vec3 b = cross(v, vec3(0, 0, 1));
    if (dot(b, b) < 0.01) b = cross(v, vec3(0, 1, 0));
    return b;
}

bool IntersectCylinder(vec3 origin, vec3 dir, out float t){
    vec3 A = gEndpoints[1]; vec3 B = gEndpoints[2];
    float Epsilon = 0.0000001;
    float extent = distance(A, B);
    vec3 W = (B - A) / extent;
    vec3 U = perp(W);
    vec3 V = cross(U, W);
    U = normalize(cross(V, W));
    V = normalize(V);
    float rSqr = Radius*Radius;
    vec3 diff = origin - 0.5 * (A + B);
    mat3 basis = mat3(U, V, W);
    vec3 P = diff * basis;
    float dz = dot(W, dir);
    if (abs(dz) >= 1.0 - Epsilon) {
        float radialSqrDist = rSqr - P.x*P.x - P.y*P.y;
        if (radialSqrDist < 0.0)
            return false;
        t = (dz > 0.0 ? -P.z : P.z) + extent * 0.5;
        return true;
    }

    vec3 D = vec3(dot(U, dir), dot(V, dir), dz);
    float a0 = P.x*P.x + P.y*P.y - rSqr;
    float a1 = P.x*D.x + P.y*D.y;
    float a2 = D.x*D.x + D.y*D.y;
    float discr = a1*a1 - a0*a2;
    if (discr < 0.0)
        return false;

    if (discr > Epsilon) {
        float root = sqrt(discr);
        float inv = 1.0/a2;
        t = (-a1 + root)*inv;
        return true;
    }

    t = -a1/a2;
    return true;
}

vec3 ComputeLight(vec3 L, vec3 N, bool specular){
    float df = max(0.0,dot(N, L));
    vec3 color = df * DiffuseMaterial;
    if (df > 0.0 && specular) {
        vec3 E = vec3(0, 0, 1);
        vec3 R = reflect(L, N);
        float sf = max(0.0,dot(R, E));
        sf = pow(sf, Shininess);
        color += sf * SpecularMaterial;
    }
    return color;
}

void main(){
    vec3 rayStart = gPosition;
    vec3 rayEnd = vec3(0);
    vec3 rayDir = normalize(rayEnd - rayStart);

    if (distance(rayStart, rayEnd) < 0.1) {discard;return;}
    float d;
    if (!IntersectCylinder(rayStart, rayDir, d)) {discard;return;}
    vec3 hitPoint = rayStart + d * rayDir;
    if (dot(hitPoint - gEndpoints[1], gEndplanes[0]) < 0.0) {discard;return;}
    if (dot(hitPoint - gEndpoints[2], gEndplanes[1]) > 0.0) {discard;return;}

    // Compute a lighting normal:
    vec3 x0 = hitPoint;
    vec3 x1 = gEndpoints[1];
    vec3 x2 = gEndpoints[2];
    float length = distance(x1, x2);
    vec3 v  = (x2 - x1) / length;
    float t = dot(x0 - x1, v);
    vec3 spinePoint = x1 + t * v;
    vec3 N  = -normalize(hitPoint - spinePoint);

    // Perform lighting and write out a new depth value:
    vec3 color = AmbientMaterial + ComputeLight(LightDirection, N, true);
    vec4 ndc = Projection * vec4(hitPoint, 1);
    gl_FragDepth = ndc.z / ndc.w;
    gl_FragColor = vec4(color, 1.0);
}
