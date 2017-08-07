#version 330 core

in float fPathCoord;

out vec4 outColor;
void main(){
    outColor = vec4(fPathCoord, 0.0, 0.0, 1.0);
}

