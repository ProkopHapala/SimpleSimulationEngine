// from : http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

#version 330 core

layout(location = 0) out vec3 color;

uniform vec2  resolution;
in vec2 UV;                 // Interpolated values from the vertex shaders
uniform sampler2D texZ;     // Values that stay constant for the whole mesh.
uniform sampler2D texRGB;   // Values that stay constant for the whole mesh.

float occlusion( vec2 uv, vec2 duv, float z ){
	float z00 = textureLod( texZ, uv+duv.xy, 0 ).r;
	float z11 = textureLod( texZ, uv-duv.xy, 0 ).r;
	duv = vec2(duv.y,-duv.x);
	float z10 = textureLod( texZ, uv+duv, 0 ).r;
	float z01 = textureLod( texZ, uv-duv, 0 ).r;
	float z1 = z00 + z11;
 	float z2 = z10 + z01;
	//return (z1+z2)*0.25;
	return ( z - min( z1, z2 )*0.5 )/dot(duv,duv);
	//return max( z1, z2 )*0.5;
}

void main(){
	
	vec2 uv = (gl_FragCoord.xy / resolution).xy;
	float z = textureLod( texZ, uv, 0 ).r; 

	//float dzsum = 0;
	//float z    =     textureLod( texZ, uv, 0 ).r;
	//float dz   = z - textureLod( texZ, uv, 3 ).r;  dzsum += dz;

	//float d = 0.02 + 0.008*sin(uv.x*uv.y*48145.0),  w  = 0.0;
	//float d = 0.02 + 0.008*sin(z*48145989.0),  w  = 0.0;
	//float d = 0.02 + 0.008*sin(z*481009.0),  w  = 0.0;
	float d = 0.02 + 0.008*sin(z*881009.0),  w  = 0.0;
	//duv = d*fract(sin(1548.0*vec2( uv.x*uv.y, 100* uv.x-uv.y ))*43758.5453);
	//float phi = 1548748.0 * uv.x*uv.y;	vec2 duv = 2*d * vec2( sin( phi ), cos( phi ) ); 
	//duv = vec2( 1.4, 1.4 );
 
	float occ;
	float occ1 = occlusion( uv, vec2( 0.5*d,  0.0   ), z )*1.0;
	float occ2 = occlusion( uv, vec2( 1.4*d,  1.4*d ), z )*2.0;
	float occ3 = occlusion( uv, vec2( 4.0*d,  0.0   ), z )*3.0*0.0;
	//occ  = max( occ1, occ2 ); 
	//occ  = max( occ , occ3 ); 
	occ = occ1 + occ2 + occ3;
	//occ = occ1;

	if( occ > 0 ){
		occ = 0.6*occ/(0.2 + occ );
		vec3 zcolor = vec3( occ, occ, occ );
		color       = textureLod( texRGB, uv, 0 ).rgb - zcolor;
		//color     =  zcolor;
	}else{
		color       = textureLod( texRGB, uv, 0 ).rgb;
	}
	
}


