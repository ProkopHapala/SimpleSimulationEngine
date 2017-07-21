#version 330 core

layout(location = 0) in vec3 model_vpos;
layout(location = 1) in vec3 model_vnor;
layout(location = 2) in vec3 pose_pos;
layout(location = 3) in vec3 pose_dir; 
layout(location = 4) in vec3 pose_Up;  
layout(location = 5) in vec3 pose_sc;  

uniform vec3 camPos;
uniform mat4 camMat;

noperspective out vec3 vnor_world;

void main(){
    mat3 rotMat  = mat3( cross(pose_dir,pose_Up), pose_Up, pose_dir );
    vec3 world_vpos = rotMat * (model_vpos*pose_sc) + pose_pos;
    //vec3 world_vpos = rotMat * model_vpos + pose_pos;
    vnor_world      = rotMat *  model_vnor;
    //vec3 world_vpos = model_vpos + pose_pos;
    //vnor_world      = model_vnor;
    gl_Position     = camMat * vec4( world_vpos-camPos, 1.0 );
}

