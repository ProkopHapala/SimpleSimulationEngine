

uniform vec2  resolution;
uniform float time;
uniform mat3  camMat;

/*
vec3 accord_1( vec2 p ){
    p *= 6.0;
    vec2 sp  = sin(p); vec2 cp  = cos(p);
    vec2 spn = sp;     vec2 cpn = cp;
    vec3 fd = vec3(
        cpn.x*spn.y,
        spn.x*cpn.y,
        spn.x*spn.y
    );
    for(int i=1; i<4; i++){
        vec2 cpn_;
        cpn_ = cpn*cp - spn*sp;
        spn  = cpn*sp + spn*cp;
        cpn=cpn_;
        fd -= vec3(  
            cpn.x*spn.y,
            spn.x*cpn.y,
            spn.x*spn.y /float(i+1)
        );
    }
    fd.xy*=6.0;
    return fd;
    
}
*/

vec2 mul_complex( vec2 a, vec2 b ){ return vec2( a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x ); }
mat2 mul_complex( mat2 a, mat2 b ){ return mat2( a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0] ); }

vec3 sin_fd( vec2 p, vec2 freq ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p) * freq;
    return vec3( cp.x*sp.y, sp.x*cp.y, sp.x*sp.y );
}

vec3 lorenz2D( vec2 sp, vec2 cp, float w2 ){
    float D = 1.0/(w2 + dot(sp,sp) );
    return vec3( -2.0*D*D*cp*sp, D )*w2;
}

vec3 sin_poles( vec2 p, vec2 freq, float w2 ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p);
    return lorenz2D( sp, cp*freq, w2 );
}

vec3 poles_harmonics( vec2 p ){
    vec2 freq = vec2(1.0); 
    float w2 = 1.5;
    p       *= freq;
    //vec2 sp  = sin(p);
    //vec2 cp  = cos(p);
    mat2 cs  = mat2( cos(p), sin(p) ); 
    mat2 csn = cs;
    vec3 fd  = lorenz2D( csn[1], csn[0]*freq, w2 ); 
    // harmonic series 1,2,3,4
    //csn = mul_complex( csn, cs );    fd += lorenz2D( csn[1], csn[0]*freq*2.0, w2 )*-0.5;
    //csn = mul_complex( csn, cs );    fd += lorenz2D( csn[1], csn[0]*freq*3.0, w2 )* 0.333;
    //csn = mul_complex( csn, cs );    fd += lorenz2D( csn[1], csn[0]*freq*4.0, w2 )*-0.25;   
    // geometric series 1,2,4,8
    csn = mul_complex( csn, csn );    fd += lorenz2D( csn[1], csn[0]*freq*2.0, w2 )*-1.25;
    csn = mul_complex( csn, csn );    fd += lorenz2D( csn[1], csn[0]*freq*4.0, w2 )* 0.75;
    csn = mul_complex( csn, csn );    fd += lorenz2D( csn[1], csn[0]*freq*8.0, w2 )*-0.55;  
    return fd*2.0;
}

/*
vec3 rot_sin( vec2 ex, vec2 ey, vec2 a ){
    // exp(i*k.p) = exp(i*(kx*px+ky.py)) = exp(i*kx*px) * exp(i*ky*py)
    vec2 exy = mul_complex( ex, ey );
    float s = dot(a,sp);
    return vec3( 0.0,0.0, s*s-0.5 );
}

vec3 sin_ndir( vec2 p ){
    vec2 freq = vec2(1.0); 
    p       *= freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p);
    //vec3 fd  = rot_sin( cp, sp, vec2(0.70710678118,-0.70710678118) );
    vec3 fd  = rot_sin( cp, sp, vec2(1.0,0.0) );
    return fd;
}
*/


vec3 func( vec2 p ){
    //return sin_fd( p, vec2(2.0) ) - sin_fd( p, vec2(4.0) )*0.2;
    //return sin_poles( p, vec2(2.0), 0.5 ) - sin_poles( p, vec2(4.0), 0.25 )*0.2;
    return poles_harmonics( p );
    //return sin_ndir( p );
}

vec4 checkFuncDeriv( vec2 p ){
    vec3 fd = func( p );
    if(p.x>0.0){
        float dd = 0.001;
        fd.xy=vec2(
            func( p + vec2(dd,0.0) ).z - fd.z ,
            func( p + vec2(0.0,dd) ).z - fd.z
        )/dd;
    }
    return vec4( vec3(fd.xy,0.0)*0.1 + vec3(fd.z*0.5+0.5), 0.0);
    //return vec4( vec3(fd.z*0.5+0.5), 0.0);
    //return vec4( vec3(fd.xy,0.0)*0.1 + 0.5, 0.0);
}

void main( ){	
    
	vec2 p    = 2.0 * gl_FragCoord.xy / resolution.xy - 1.0;
	p *=6.0;	
    gl_FragColor = checkFuncDeriv(p); 
}
