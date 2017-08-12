//==========================================
//>>const3D.vert
#version 330 core
layout(location = 0) in vec3 vertPos_model;
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;
void main(){
	vec3 position_world = modelPos + modelMat * vertPos_model;
	gl_Position         = camMat   * vec4( position_world-camPos, 1 );
}
//<<
//==========================================
//>>const3D.frag
#version 330 core
smooth in vec4 gl_FragCoord;
out vec3 color;
uniform vec4 baseColor;
void main(){
	color = baseColor.xyz;
	gl_FragDepth = gl_FragCoord.z;
}
//<<
//>>color3D.frag
#version 330 core
smooth in vec4 gl_FragCoord;
smooth in vec3 fragColor;
out vec3 color;
uniform vec3 light_dir;
void main(){
	color = fragColor;
	gl_FragDepth = gl_FragCoord.z;
}﻿
//<<
//==========================================
//>>color3D.vert
#version 330 core
layout(location = 0) in vec3 vertPos_model;
layout(location = 1) in vec3 vertColor;
smooth out vec3 fragColor;
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;
void main(){
	vec3 position_world = modelPos + modelMat * vertPos_model;
	gl_Position         = camMat   * vec4( position_world-camPos, 1 );
	fragColor           = vertColor;
}
//<<
//==========================================
//>>shade3D.vert
#version 330 core
layout(location = 0) in vec3 vertPos_model;
layout(location = 1) in vec3 vertNormal_model;
noperspective out vec3 fragNormal_world;
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;
void main(){
	vec3 position_world = modelPos + modelMat * vertPos_model;
	gl_Position         = camMat   * vec4( position_world-camPos, 1 );
	fragNormal_world    = modelMat * vertNormal_model;
}
//<<
//==========================================
//>>shade3D.frag
#version 330 core
smooth        in vec4 gl_FragCoord;
noperspective in vec3 fragNormal_world;
out vec3 color;
uniform vec3 cam_pos;   // world - but does not realy matter
uniform vec3 light_pos; // world - but does not realy matter
uniform	vec3 lightColor;
uniform	vec3 diffuseColor;
uniform	vec3 ambientColor;
uniform	vec3 specularColor;
void main(){
	// difuse
	vec3 n = normalize( fragNormal_world );
	vec3 l = normalize( -light_pos );
	float cosTheta = clamp( dot( n,l ), 0,1 );
	// specular
	vec3 E = normalize(-cam_pos);
	vec3 R = reflect(-l,n);
	float cosAlpha = clamp( dot( E,R ), 0,1 );
	color = ambientColor + lightColor*( diffuseColor*cosTheta  +  specularColor*pow(cosAlpha,5) );
}
//<<
