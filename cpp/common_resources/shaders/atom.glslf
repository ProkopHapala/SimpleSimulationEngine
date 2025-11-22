varying vec3 vColor;
varying vec2 vUv;
varying vec3 vViewPosition;
varying float vRadius;

void main() {
    // 1. Calculate distance from center (0.5, 0.5)
    vec2 uvOffset = vUv - 0.5;
    float distSq = dot(uvOffset, uvOffset);
    
    // 2. Discard if outside circle (radius 0.5)
    if (distSq > 0.25) discard;

    // 3. Calculate normal for lighting (Impostor)
    // We need to reconstruct the Z component of the sphere surface
    // x^2 + y^2 + z^2 = r^2  =>  z = sqrt(r^2 - x^2 - y^2)
    // Here r=0.5. uvOffset is x,y.
    
    float z = sqrt(0.25 - distSq);
    vec3 normal = normalize(vec3(uvOffset.x, uvOffset.y, z));

    // 4. Simple Lighting (Directional)
    vec3 lightDir = normalize(vec3(0.5, 0.5, 1.0));
    float diffuse = max(dot(normal, lightDir), 0.0);
    
    // Ambient
    vec3 ambient = vec3(0.3);
    
    // Specular (Phong)
    vec3 viewDir = normalize(vec3(0.0, 0.0, 1.0)); // In tangent space/billboard space, view is roughly +Z
    vec3 reflectDir = reflect(-lightDir, normal);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32.0);
    vec3 specular = vec3(0.5) * spec;

    // Combine
    vec3 finalColor = vColor * (ambient + diffuse) + specular;

    gl_FragColor = vec4(finalColor, 1.0);
}
