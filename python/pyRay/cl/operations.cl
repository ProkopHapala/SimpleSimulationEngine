
float opU( float d1, float d2 ){
	return fmin(d1,d2);
}

float opS( float d1, float d2 ){
    return fmax(-d2,d1);
}

float opI( float d1, float d2 ){
    return fmax(d1,d2);
}

float3 opRep( float3 p, float3 c ){
    return fmod( p, c)-0.5f*c;
}

float3 opRepRot( float3 p, float times ){
    float phi  = atan2(p.z,p.x);
    float dphi = 3.14159265359f / times;
    phi        = fmod( phi+dphi, 2.0f*dphi )-dphi;  
    float r    = length(p.xz);
    return     (float3)( r*cos(phi), p.y, r*sin(phi) );
}

/*
float3 opTwist( float3 p ){
    float  c = cos(10.0f*p.y+10.0f);
    float  s = sin(10.0f*p.y+10.0f);
    mat2   m = mat2(c,-s,s,c);
    return (float3)(m*p.xz,p.y);
}
*/
