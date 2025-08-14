#version 330 core
in vec2 uv;
out vec4 fragColor;
void main(){
    vec2 p = uv - vec2(0.5);
    float d2 = dot(p,p);
    float a = exp(-40.0*d2); // Gaussian blob
    fragColor = vec4(vec3(a), 1.0);
}
