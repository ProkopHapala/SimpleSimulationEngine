
uniform vec2  resolution;
uniform float time;
uniform mat3  camMat;

#define MAX_BOUNCES 2
float gamma = 2.2;

struct DirectionalLight{  vec3 d; vec3 c; };
struct Material{ vec3 color; float gloss; };

DirectionalLight sunLight = DirectionalLight( normalize(vec3(1.0, 0.5, 0.5)), vec3(1.0) );

float Lorenz( float x ){ return 1.0/(1.0+x*x); }
    
void main( ){	
    mat3 camMat_ = camMat;
	vec3 uvz = vec3( 2.0 * gl_FragCoord.xy / resolution.xy - 1.0, 5.0 );
	
	vec3 p  = vec3(0.0, 0.0, 10.0 );
    vec3 d  = normalize(vec3(resolution.x/resolution.y * uvz.x, uvz.y, -uvz.z ) );
    Ray ray = Ray(camMat_*p, camMat_*d);
 
	vec4 hit = scene(ray);

	Material mat = Material( vec3(1.0,0.5,0.5), 1.0 ); // we don't care for now

    if( hit.x<SKY_DIST ){
        float c_diffuse  = clamp( dot(hit.yzw,sunLight.d), 0.0, 1.0);
        //-- specular
        vec3  nn = sunLight.d-ray.d;
        float cn = dot(hit.yzw,nn);
        float c_specular = Lorenz( 100.0*(1.0-clamp( (cn*cn)/dot(nn,nn), 0.0, 1.0)) );
        //-- output
        //gl_FragColor = vec4( (c_diffuse + c_specular*mat.gloss)*mat.color + vec3(0.1,0.1,0.2)*mat.color, 1.0 ); 
        //gl_FragColor = vec4( (hit.yzw+vec3(1.0))*0.5, 1.0 );   
        //gl_FragColor = vec4( vec3(log(hit.x)-2.0), 1.0 );   
        //gl_FragColor = vec4( sin(hit.x*4.0),sin(hit.x),sin(hit.x*16.0), 1.0 );
        gl_FragColor = OUTPUT_PIXEL 
    }else{
        discard;
    }
    
}
