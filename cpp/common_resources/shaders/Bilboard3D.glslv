#version 330 core

layout(location = 0) in vec2 model_uv;
layout(location = 1) in vec3 pose_pos;
layout(location = 2) in vec2 pose_sc;  

uniform vec3 camPos;
uniform mat4 camMat;
uniform mat3 camRot;

noperspective out vec3 fpos_world;
//noperspective out vec2 fUV;
smooth out vec2 fUV;
noperspective out vec4 obj;

void main(){
    vec2 p2d          = model_uv*pose_sc;
    vec3 world_vpos   = pose_pos + camRot[0]*p2d.x + camRot[1]*p2d.y;

    fUV               = model_uv;
    fpos_world        = world_vpos;
    obj               = vec4( pose_pos, pose_sc.x );
    gl_Position       = camMat * vec4( world_vpos-camPos, 1.0 );
    
    //gl_Position = vec4( model_uv, 0.5, -0.5  );
}

