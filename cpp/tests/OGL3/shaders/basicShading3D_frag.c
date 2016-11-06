#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/
// http://www.geeks3d.com/20130514/opengl-interpolation-qualifiers-glsl-tutorial/
// interpolation qualifiers : flat, smooth, noperspective

// IN --- Interpolated values from the vertex shaders
smooth        in vec4 gl_FragCoord;
noperspective in vec3 fragNormal_world;

// OUT --- Ouput data
out vec3 color;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 cam_pos;   // world - but does not realy matter
uniform vec3 light_pos; // world - but does not realy matter

uniform	vec3 lightColor;
uniform	vec3 diffuseColor;
uniform	vec3 ambientColor;
uniform	vec3 specularColor;

void main(){

	//gl_FragDepth = gl_FragCoord.z;

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
