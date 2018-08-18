#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in        vec4 gl_FragCoord;
noperspective in vec3 world_pos;
noperspective in vec3 world_nor;
//noperspective in vec4 tx_;

// OUT --- Ouput data
out vec3 color;

uniform vec3 camPos;
// UNI --- Values that stay constant for the whole mesh.
//uniform vec3 camPos;   // world - but does not realy matter
uniform vec3 lightPos; // world - but does not realy matter

uniform	vec3 lightColor;
uniform	vec3 diffuseColor;
uniform	vec3 ambientColor;
uniform	vec3 specularColor;

void main(){

	float iso = step( 0.9, fract(world_pos.y) ); 
	float h = world_pos.y*0.05;
	//vec3 d=world_pos-camPos;
	//float l  = ( log( dot(d,d) )-6.0 )*0.4;
	//color = vec3(h,0.5-0.5*h,iso);
	//color = vec3(l*c);

    //color = sin(world_nor*10000.0);
    
    //color = vec3( sin(tx_.a*100.0) );
    //color = sin(world_pos*10.0);
    //color = vec3( sin(world_pos.y*10.0) );
    //color = vec3( sin( tx_.x*15.0 )  );
    //color = vec3( sin( tx_.rgb*150.0 )  );
    //color = vec3( sin( world_nor*150.0 )  );
    //color = world_nor;

    // difuse
	vec3 l = normalize( world_pos-lightPos );
	vec3 E = normalize( camPos-world_pos   );
	vec3 n = normalize( world_nor );
    if( dot(E,n)<0 ) n*=-1.0;
	float cosTheta = clamp( -dot( n,l ), 0.0,1.0 );
	
	// specular 
	vec3 R = reflect( l,n);
	float cosAlpha = clamp( dot( E,R ), 0.0,1.0 );
	
	color = ambientColor*(1+0.5*fract(world_pos.y*1.0) ) + lightColor*( diffuseColor*cosTheta  +  specularColor*pow(cosAlpha,16.0) );
	
	//color = sin(world_pos*1.0);
	//color = vec3( sin(world_pos.y*1.0) , log( world_pos.y)-2,  log( world_pos.y)-3 );
	
	//color    = world_nor * 10.0;  color.y  = 0.0;  color   += vec3(0.5);
	
	//color = world_nor;
	
	//gl_FragDepth = -dist;
	//gl_FragDepth = gl_FragCoord.z;
}
