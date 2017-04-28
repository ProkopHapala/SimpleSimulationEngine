
#define PI 3.14159275358979
#define POSITIVE_INF 1e+8
#define SKY_DIST     1e+7

// Structs
struct Ray     { vec3 o;  vec3 d; };
struct Plane   { vec3  n;         float C; };
struct Slab    { vec3  n;         float Cmin; float Cmax; };
struct Sphere  { vec3  p;         float R; };
struct Cylinder{ vec3  p; vec3 d; float R; };
struct Cone    { vec3  c; float cosa; vec3  v; float h; };

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

// === Cone

#define ADD( SURF ) { ts1 = ray_ts( SURF, ray ); if( ts1.x < hit.x ){ hit = vec4(ts1.x, normal(SURF,point(ray,ts1.x))); } }
#define SUB( SURF ) { vec2 ts2 = ray_ts( SURF, ray ); if( (ts1.x>ts2.x) && (ts1.x<ts2.y) ){ if( ts2.y<ts1.y ){ hit = vec4( ts2.y, normal(SURF,point(ray,ts2.y))*-1.0 ); }else{ hit = vec4( POSITIVE_INF, vec3(0.0) ); } } }
