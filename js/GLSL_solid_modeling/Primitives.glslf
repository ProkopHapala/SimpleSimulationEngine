
#define PI 3.14159275358979
#define POSITIVE_INF 1e+8
#define SKY_DIST     1e+7
#define PREC_STEP    1e-5

// Structs
struct Ray     { vec3 o;  vec3 d; };
struct Plane   { vec3  n;          float C; };
struct Slab    { vec3  n;          float Cmin; float Cmax; };
struct Sphere  { vec3  p;          float R; };
struct Tube    { vec3  p; vec3 d;  float R; };                 // infinite cylinder
struct Cylinder{ vec3  a; vec3 b;  float R; };                 // finite cylinder
struct Cone    { vec3  p; vec3 d;  float h; float cosa; };

// === Ray 

vec3 point( Ray r, float t ){ return r.o + t*r.d; }

// === Plane

float dist( Plane pl, vec3 p ){
    float c = dot(pl.n,p);
    return c - pl.C;
}

float ray_t(Plane pl, Ray ray){
    float cnd = dot(pl.n, ray.d);
    float c   = dot(pl.n, ray.o);
    return -(c + pl.C) / cnd;
}

// === Slab

float dist( Slab sl, vec3 p ){
    float c = dot(sl.n,p);
    if (c<sl.Cmin){ return sl.Cmin-c; } else if (c>sl.Cmax){ return c-sl.Cmax; } else { return -1.0; };
}

vec2 ray_ts( Slab sl, Ray ray){
    float cnd = dot(sl.n, ray.d);
    float c   = dot(sl.n, ray.o);
    float t1  = -(c + sl.Cmin) / cnd;
    float t2  = -(c + sl.Cmax) / cnd;
    if(t1<t2){ return vec2(t1,t2); }else{ return vec2(t2,t1); };
}

float ray_t( Slab sl, Ray ray){
    vec2 ts = ray_ts( sl, ray);
    if (ts.x>0.0){ return ts.x; }else{ return ts.y; };
}

vec3 normal( Slab sl, vec3 p ){ 
    float c = dot(sl.n,p);
    if( c > 0.5*(sl.Cmin+sl.Cmax) ){ return sl.n; }else{ return -sl.n; } 
}

// === Sphere

float dist2( Sphere sph, vec3 p ){
    vec3 dp = p - sph.p;
    return dot(dp,dp)-(sph.R*sph.R);
}

float ray_t( Sphere sph, Ray ray ){
    vec3 op   = sph.p - ray.o;
    float b   = dot(op, ray.d);
    float det = b*b - dot(op,op) + sph.R*sph.R;
    if (det<0.0) return POSITIVE_INF;
    det       = sqrt(det);
    float t   = b - det; 
    if (t < 0.0) t = b + det;
    return t;
}

vec2 ray_ts( Sphere sph, Ray ray ){
    vec3 op   = sph.p - ray.o;
    float b   = dot(op, ray.d);
    float det = b*b - dot(op,op) + sph.R*sph.R;
    if (det<0.0) return vec2(POSITIVE_INF,POSITIVE_INF);
    det       = sqrt(det);
    return vec2( b-det, b+det );
}

vec3 normal( Sphere sph, vec3 p ){ return (p-sph.p)/sph.R; }

// === Cylinder

float dist2( Cylinder cl, vec3 p ){
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    vec3  da  = p-cl.a;
    vec3  db  = p-cl.b;
    vec3  aXb = cross(da,db);
    vec3  dc  = cl.b-cl.a;
    float d2  = dot(aXb,aXb)/dot( dc,dc );
    return sqrt(d2) - cl.R;
}

vec2 ray_ts( Cylinder cl, Ray ray ){

    vec3  cc  = 0.5*(cl.a+cl.b);
    float ch  = length(cl.b-cl.a);
    vec3  ca  = (cl.b-cl.a)/ch;
    ch       *= 0.5;

    vec3  oc = ray.o - cc;

    float card = dot(ca,ray.d);
    float caoc = dot(ca,oc);
    
    float a = 1.0 - card*card;
    float b = dot( oc, ray.d) - caoc*card;
    float c = dot( oc, oc   ) - caoc*caoc - cl.R*cl.R;
    float h = b*b - a*c;
    if( h<0.0 ) return vec2(POSITIVE_INF,POSITIVE_INF);
    h = sqrt(h);
    
    vec2 ts;
    
    float t = (-b-h)/a;
    float y = caoc + t*card; // position of hit along axis
    if( abs(y)<ch ){
        ts.x = t;
    } else{
        float t = (sign(y)*ch - caoc)/card;
        if( abs(b+a*t)<h ){
            ts.x = t;
        }else{
            ts.x = POSITIVE_INF;
        }
    }
    
    if( ts.x < POSITIVE_INF ){
        float t = (-b+h)/a;
        float y = caoc + t*card; // position of hit along axis
        if( abs(y)<ch ){
            ts.y = t;
        } else{
            ts.y = (sign(y)*ch - caoc)/card;
        }
    }
    
    return ts;
}

vec3 normal( Cylinder cl, vec3 p ){ 
    vec3  dc  = cl.b-cl.a;
    float h   = length(dc); 
    dc       *= (1.0/h);
    vec3  dp  = p-cl.a;
    float c   = dot(dp,dc);
    if      (c<PREC_STEP){
        return  dc*-1.0;
    }else if(c>h-PREC_STEP){
        return  dc;
    }else {
        return normalize( dp - dc*c );
    }
}

// === Cone

vec2 ray_ts(Cone cone, Ray ray){
    
    /*
    float c1  = dot(ray.d,cone.d);
    float cs2 = cone.cosa*cone.cosa;
    
    float c   = dot(ray.o-cone.p,cone.d);
    float t   = -(c - cone.h)/c1;
    //float t   = -(c + cone.h)/c1;
    vec3 p    = ray.o + t*ray.d;
    vec3 d    = p - (cone.p+(cone.d*cone.h));
    if( dot(d,d) < ((1.0-cs2)/cs2)*cone.h*cone.h ){
        return vec2(t,t);
    }else{
        return vec2(POSITIVE_INF,POSITIVE_INF);
    }
    */
    
    vec3  co  = ray.o -   cone.p;
    float c1  = dot(ray.d,cone.d);
    float c2  = dot(co   ,cone.d);
    float cs2 = cone.cosa*cone.cosa;
    float a   = c1*c1 - cs2;
    float b   = c1*c2 - cs2*dot(ray.d,co);
    float c   = c2*c2 - cs2*dot(co,co);
    float det = b*b - a*c;
    if (det < 0.) return vec2(POSITIVE_INF,POSITIVE_INF);
    det       = sqrt(det);
    float t1  = (-b - det) / a;
    float t2  = (-b + det) / a;
    vec2 ts = vec2(t1,t2);
    
    vec3  p = ray.o + t1*ray.d;
    float h = dot(p - cone.p, cone.d);
    
    if(h<0.0)    return vec2(POSITIVE_INF,POSITIVE_INF);
    if(h>cone.h){
        float c   = dot(ray.o-cone.p,cone.d);
        float t   = -(c - cone.h)/c1;
        vec3 p    = ray.o + t*ray.d;
        vec3 d    = p - (cone.p+(cone.d*cone.h));
        if( dot(d,d) < ((1.0-cs2)/cs2)*cone.h*cone.h ){
            return vec2(t,t);
        }else{
            return vec2(POSITIVE_INF,POSITIVE_INF);
        }
    }
    return ts;
  
}

vec3 normal( Cone cone, vec3 p ){ 
    vec3   dp = p-cone.p;
    if( dot(cone.d,dp) > (cone.h-PREC_STEP) ){
        //return -cone.d;
        return dp;
    }else{
        float cp = dot(cone.d,dp)/dot(dp,dp);
        return normalize( cone.d - cp*dp );
    }
}

#define ADD( SURF ) { ts1 = ray_ts( SURF, ray ); if( ts1.x < hit.x ){ hit = vec4(ts1.x, normal(SURF,point(ray,ts1.x))); } }
#define SUB( SURF ) { vec2 ts2 = ray_ts( SURF, ray ); if( (ts1.x>ts2.x) && (ts1.x<ts2.y) ){ if( ts2.y<ts1.y ){ hit = vec4( ts2.y, normal(SURF,point(ray,ts2.y))*-1.0 ); }else{ hit = vec4( POSITIVE_INF, vec3(0.0) ); } } }
