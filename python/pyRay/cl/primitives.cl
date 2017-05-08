float sdPlaneY( float3 p ){
	return p.y;
}

float sdPlane( float3 p, float3 dir, float C ){
	return dot(dir,p) - C;
}

float sdSphere( float3 p, float s ){
    return length(p)-s;
}

float sdBox( float3 p, float3 b ){
    float3 d = fabs(p) - b;
    return fmin(fmax(d.x,fmax(d.y,d.z)),0.0f) + length(fmax(d,0.0f));
}

float sdEllipsoid( float3 p, float3 r ){
    return (length( p/r ) - 1.0f) * fmin(fmin(r.x,r.y),r.z);
}

float udRoundBox( float3 p, float3 b, float r ){
    return length(fmax(fabs(p)-b,0.0f))-r;
}

float sdTorus( float3 p, float2 t ){
    return length( (float2)(length(p.xz)-t.x,p.y) )-t.y;
}

float sdHexPrism( float3 p, float2 h ){
    float3 q = fabs(p);
    float d1 = q.z-h.y;
    float d2 = fmax((q.x*0.866025f+q.y*0.5f),q.y)-h.x;
    return length(fmax((float2)(d1,d2),0.0f)) + fmin(fmax(d1,d2), 0.0f );
}

float sdCapsule( float3 p, float3 a, float3 b, float r ){
	float3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0f, 1.0f );
	return length( pa - ba*h ) - r;
}

float sdTriPrism( float3 p, float2 h ){
    float3 q = fabs(p);
    float d1 = q.z-h.y;
    float d2 = fmax(q.x*0.866025f+p.y*0.5f,-p.y)-h.x*0.5f;
    return length(fmax((float2)(d1,d2),0.0f)) + fmin(fmax(d1,d2), 0.0f );
}

float sdCylinder( float3 p, float2 h ){
  float2 d = fabs((float2)(length(p.xz),p.y)) - h;
  return fmin(fmax(d.x,d.y),0.0f) + length(fmax(d,0.0f));
}

float sdOgiv( float3 p, float2 h ){
  float R = (h.x*h.x + h.y*h.y)/(2*h.y);
  return  length(p.xz) - sqrt( R*R - p.y*p.y ) + (R-h.y);
}

float sdCone( float3 p, float3 c ){
    float2 q = (float2)( length(p.xz), p.y );
    float d1 = -q.y-c.z;
    float d2 = fmax( dot(q,c.xy), q.y);
    return length(fmax((float2)(d1,d2),0.0f)) + fmin(fmax(d1,d2), 0.0f );
}

float sdConeSection( float3 p, float h, float r1, float r2 ){
    float d1 = -p.y - h;
    float q = p.y - h;
    float si = 0.5f*(r1-r2)/h;
    float d2 = fmax( sqrt( dot(p.xz,p.xz)*(1.0f-si*si)) + q*si - r2, q );
    return length(fmax((float2)(d1,d2),0.0f)) + fmin(fmax(d1,d2), 0.0f );
}

float sdPryamid4(float3 p, float3 h ) {     // h = { cos a, sin a, height }
    // Tetrahedron = Octahedron - Cube
    float box = sdBox( p - (float3)(0.0f,-2.0f*h.z,0.0f), (float3)(2.0f*h.z) );
    float d = 0.0f;
    d = fmax( d, fabs( dot(p, (float3)( -h.x, h.y, 0.0f )) ));
    d = fmax( d, fabs( dot(p, (float3)(  h.x, h.y, 0.0f )) ));
    d = fmax( d, fabs( dot(p, (float3)(  0.0f, h.y, h.x )) ));
    d = fmax( d, fabs( dot(p, (float3)(  0.0f, h.y,-h.x )) ));
    float octa = d - h.z;
    return fmax(-box,octa); // Subtraction
 }

float length2( float2 p ){
	return sqrt( p.x*p.x + p.y*p.y );
}

float length6( float2 p ){
	p = p*p*p; p = p*p;
	return pow( p.x + p.y, 1.0f/6.0f );
}

float length8( float2 p ){
	p = p*p; p = p*p; p = p*p;
	return pow( p.x + p.y, 1.0f/8.0f );
}

float sdTorus82( float3 p, float2 t ){
    float2 q = (float2)(length2(p.xz)-t.x,p.y);
    return length8(q)-t.y;
}

float sdTorus88( float3 p, float2 t ){
    float2 q = (float2)(length8(p.xz)-t.x,p.y);
    return length8(q)-t.y;
}

float sdCylinder6( float3 p, float2 h ){
    return fmax( length6(p.xz)-h.x, fabs(p.y)-h.y );
}

