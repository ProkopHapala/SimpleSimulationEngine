#version 330 core

flat in vec3  p_start;
flat in vec3  p_end;
in vec3  fpos_world;
in float fPathCoord;

out vec4 outColor;

uniform vec3 camPos;

vec4 iCylinder( in vec3 ro, in vec3 rd,  in vec3 pa, in vec3 pb, float ra ){
    // from https://www.shadertoy.com/view/4lcSRn
    vec3  cc   = 0.5*(pa+pb);
    float ch   = length(pb-pa);
    vec3  ca   = (pb-pa)/ch;
    ch        *= 0.5;
    vec3  oc   = ro - cc;
    float card = dot(ca,rd);
    float caoc = dot(ca,oc);
    float a    = 1.0 - card*card;
    float b    = dot( oc, rd) - caoc*card;
    float c    = dot( oc, oc) - caoc*caoc - ra*ra;
    float h    = b*b - a*c;
    if( h<0.0 ) return vec4(-1.0);
    h = sqrt(h);
    float t1 = (-b-h)/a;
    //float t2 = (-b+h)/a;
    float y = caoc + t1*card;
    // body
    if( abs(y)<ch ) return vec4( t1, normalize( oc+t1*rd - ca*y ) );
    // caps
    float sy = sign(y);
    float tp = (sy*ch - caoc)/card;
    if( abs(b+a*tp)<h ){ return vec4( tp, ca*sy ); }
    return vec4(-1.0);
}


void main(){

    vec3 ray0 = camPos; 
    vec3 hRay = fpos_world - camPos;
    hRay = normalize(hRay);

    // raytrace
	vec4 tnor = iCylinder( ray0, hRay, p_start, p_end, 0.10 );
	float t = tnor.x;
    
    // shading/lighting	
	vec3 col = vec3(0.0);
	if( t>0.0 ){
	    vec3  pos = ray0 + t*hRay;
		vec3  nor = tnor.yzw;
		float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
		float amb = 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0));
		col = vec3(0.2,0.3,0.4)*amb + vec3(0.8,0.7,0.5)*dif;
		
	    col = sqrt( col );
	    outColor = vec4( col, 1.0 );
	}else{
	    discard;
	}

	//outColor = vec4( sin((fpos_world)*100.0), 1.0);
	//outColor = vec4( sin((fpos_world-camPos)*100.0), 1.0);
	
    //outColor = vec4(fPathCoord, 0.0, 0.0, 1.0);
    
    gl_FragDepth = log(t)*0.1;
}

