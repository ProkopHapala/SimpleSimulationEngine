


// More info here: http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
// http://iquilezles.org/www/articles/raymarchingdf/raymarchingdf.htm


#define AA 1   // make this 1 is your machine is too slow

//------------------------------------------------------------------

float sdPlane    ( vec3 p )                  { return p.y; }
float sdSphere   ( vec3 p, float s )         { return length(p)-s; }
float sdBox      ( vec3 p, vec3 b )          { vec3 d = abs(p) - b; return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0)); }
float sdEllipsoid( in vec3 p, in vec3 r )    { return (length( p/r ) - 1.0) * min(min(r.x,r.y),r.z); }
float udRoundBox ( vec3 p, vec3 b, float r ) { return length(max(abs(p)-b,0.0))-r; }
float sdTorus    ( vec3 p, vec2 t )          { return length( vec2(length(p.xz)-t.x,p.y) )-t.y; }
float length2    ( vec2 p )                  { return sqrt( p.x*p.x + p.y*p.y ); }
float length6    ( vec2 p )                  { p = p*p*p; p = p*p; return pow( p.x + p.y, 1.0/6.0 ); }
float length8    ( vec2 p )                  { p = p*p; p = p*p; p = p*p; return pow( p.x + p.y, 1.0/8.0 ); }
float sdTorus82  ( vec3 p, vec2 t )          { vec2 q = vec2(length2(p.xz)-t.x,p.y); return length8(q)-t.y; }
float sdTorus88  ( vec3 p, vec2 t )          { vec2 q = vec2(length8(p.xz)-t.x,p.y); return length8(q)-t.y; }
float sdCylinder6( vec3 p, vec2 h )          { return max( length6(p.xz)-h.x, abs(p.y)-h.y ); }

float sdHexPrism( vec3 p, vec2 h ){
    vec3 q = abs(p);
#if 0
    return max(q.z-h.y,max((q.x*0.866025+q.y*0.5),q.y)-h.x);
#else
    float d1 = q.z-h.y;
    float d2 = max((q.x*0.866025+q.y*0.5),q.y)-h.x;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
#endif
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r ){
	vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h ) - r;
}

float sdTriPrism( vec3 p, vec2 h ){
    vec3 q = abs(p);
#if 0
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
#else
    float d1 = q.z-h.y;
    float d2 = max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
#endif
}

float sdCylinder( vec3 p, vec2 h ){
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdCone( in vec3 p, in vec3 c ){
    vec2 q = vec2( length(p.xz), p.y );
    float d1 = -q.y-c.z;
    float d2 = max( dot(q,c.xy), q.y);
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float sdConeSection( in vec3 p, in float h, in float r1, in float r2 ){
    float d1 = -p.y - h;
    float q = p.y - h;
    float si = 0.5*(r1-r2)/h;
    float d2 = max( sqrt( dot(p.xz,p.xz)*(1.0-si*si)) + q*si - r2, q );
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float sdPryamid4(vec3 p, vec3 h ){ 
    // h = { cos a, sin a, height }{
    // Tetrahedron = Octahedron - Cube
    float box = sdBox( p - vec3(0,-2.0*h.z,0), vec3(2.0*h.z) );
    float d = 0.0;
    d = max( d, abs( dot(p, vec3( -h.x, h.y, 0 )) ));
    d = max( d, abs( dot(p, vec3(  h.x, h.y, 0 )) ));
    d = max( d, abs( dot(p, vec3(  0, h.y, h.x )) ));
    d = max( d, abs( dot(p, vec3(  0, h.y,-h.x )) ));
    float octa = d - h.z;
    return max(-box,octa); // Subtraction
 }

//------------------------------------------------------------------

float opU( float d1, float d2 ){ return min(d1,d2);  }
float opS( float d1, float d2 ){ return max(-d2,d1); }
float opI( float d1, float d2 ){ return max(d1,d2);  }

vec3  opRep( vec3 p, vec3 c     ){ return mod(p,c)-0.5*c;        }
vec3 opTwist( vec3 p ){
    float  c = cos(10.0*p.y+10.0);
    float  s = sin(10.0*p.y+10.0);
    mat2   m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}

vec2  putObj( vec2 d1, vec2 d2  ){ return (d1.x<d2.x) ? d1 : d2; }  // combines object including material

