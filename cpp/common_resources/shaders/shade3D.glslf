#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/
// http://www.geeks3d.com/20130514/opengl-interpolation-qualifiers-glsl-tutorial/
// interpolation qualifiers : flat, smooth, noperspective

// IN --- Interpolated values from the vertex shaders
smooth        in vec4 gl_FragCoord;
//noperspective in vec3 fragNormal_world;
//noperspective in vec3 fragPos_world;
smooth in vec3 fragNormal_world;
smooth in vec3 fragPos_world;

// OUT --- Ouput data
out vec3 color;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 camPos;   // world - but does not realy matter
uniform vec3 lightPos; // world - but does not realy matter

uniform	vec3 lightColor;
uniform	vec3 diffuseColor;
uniform	vec3 ambientColor;
uniform	vec3 specularColor;

void main(){
	//gl_FragDepth = gl_FragCoord.z;

	// difuse
	vec3 l = normalize( fragPos_world-lightPos );
	vec3 E = normalize( camPos-fragPos_world   );
	vec3 n = normalize( fragNormal_world );
    if( dot(E,n)<0 ) n*=-1.0;
	float cosTheta = clamp( -dot( n,l ), 0.0,1.0 );
	
	// specular 
	vec3 R = reflect( l,n);
	float cosAlpha = clamp( dot( E,R ), 0.0,1.0 );
	
	color = ambientColor + lightColor*( diffuseColor*cosTheta  +  specularColor*pow(cosAlpha,8.0) );
	//color = ambientColor + lightColor*( diffuseColor*cosTheta );
	//color = ambientColor;
}
