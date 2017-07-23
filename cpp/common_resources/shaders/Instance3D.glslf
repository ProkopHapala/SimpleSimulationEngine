#version 330 core

//smooth        in vec4 gl_FragCoord;
noperspective in vec3 fnor_world;

out vec3 color;

void main(){
    color = fnor_world*0.5 + vec3(0.5,0.5,0.5);
}
