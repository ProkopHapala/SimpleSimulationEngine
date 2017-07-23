#version 330 core

layout(location = 0) in vec3 model_vpos;
layout(location = 1) in vec3 model_vnor;
layout(location = 2) in vec3 pose_pos;
layout(location = 3) in vec3 pose_dir; 
layout(location = 4) in vec3 pose_Up;  
layout(location = 5) in vec3 pose_sc;  

uniform vec3 camPos;
uniform mat4 camMat;

noperspective out vec3 fpos_world;
noperspective out vec3 fnor_world;
//smooth        out vec3 vpos_world;
noperspective out vec4 obj;

void main(){
    mat3 rotMat  = mat3( cross(pose_dir,pose_Up), pose_Up, pose_dir );
    vec3 world_vpos = rotMat * (model_vpos*pose_sc) + pose_pos;
    //vec3 world_vpos = rotMat * model_vpos + pose_pos;
    fnor_world        = rotMat *  model_vnor;
    fpos_world        = world_vpos;
    
    obj               = vec4( pose_pos, length(pose_sc) );
    //vec3 world_vpos = model_vpos + pose_pos;
    //vnor_world      = model_vnor;
    gl_Position     = camMat * vec4( world_vpos-camPos, 1.0 );
}

