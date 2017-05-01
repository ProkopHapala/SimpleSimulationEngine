

uniform vec2  resolution;
uniform float time;
uniform mat3  camMat;

vec3 lorenz( vec2 p, float w2 ){
	float D  = 1.0/(w2+dot(p,p));
	return vec3( p*(-2.0*D*D), D );
}

vec3 sin_fd( vec2 p, vec2 freq ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p) * freq;
    return vec3( cp.x*sp.y, sp.x*cp.y, sp.x*sp.y );
}

vec3 sin_poles( vec2 p, vec2 freq, float w2 ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p);
    float D = 1.0/(w2 + dot(sp,sp) );
    return vec3( -2.0*freq*D*D*cp*sp, D )*w2;
}

vec3 sin_poles_signed( vec2 p, vec2 freq, float w2 ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p);
    float s = cp.x*cp.y;
    float D = 1.0/(w2 + dot(sp,sp) );
    //return vec3( -2.0*freq*D*D*cp*sp*s - freq*sp*cp.yx*D, D*s )*w2;
    return vec3( (-2.0*D*cp*s - cp.yx)*freq*sp, s )*D*w2;
}

vec3 sin_net( vec2 p, vec2 freq, float w2 ){
    p *=freq;
    vec2 sp = sin(p);
    vec2 cp = cos(p);
    // --- type 1 - with singularities
    //float D  = 1.0/(w2 + sp.x*sp.y );
    //cp*=-2.0*freq*D*D;
    //return vec3( cp.x*sp.y, sp.x*cp.y, D )*w2;
    
    // --- type 2 - without nodes
    //vec2  D  = 1.0/(w2 + sp );
    //float D2 = dot(D,D);
    //cp*=-3.0*freq*D*D*D;
    //return vec3( cp.x, cp.y, D2 )*w2*w2;
    
    // --- type 3 - with nodes
    vec2  D  = 1.0/(w2 + sp*sp );
    float D2 = D.x*D.y;
    cp*=-freq*freq*D;
    return vec3( cp*sp, 1.0 )*D2*w2*w2;

}

vec3 func( vec2 p ){
/*
    return   0.5 *lorenz( p+vec2( 1.0, 0.5), 1.0 ) 
           - 1.5 *lorenz( p+vec2(-1.0,-0.5), 2.0 )
           - 0.5 *lorenz( p+vec2(-1.0,4.5), 0.5 );    
*/
//    return sin_fd( p, vec2(1.5) )*0.3;
//    return lorenz   ( p, 1.0 );
//    return sin_poles( p, vec2(1.5), 0.4 )*0.1;
//    return sin_net( p, vec2(1.5), 0.4 )*0.1;
//    return sin_poles_signed( p, vec2(1.5), 0.4 )*0.1;

//    return sin_poles_signed( p, vec2(1.5), 0.8 )*0.1 + lorenz( p, 8.0 )*8.0;
        return sin_poles_signed( p, vec2(1.0), 0.8 )*0.1 + sin_poles_signed( p, vec2(2.0), 0.4 )*0.02;
}

void move( inout vec2 p, float dt ){
    vec3 fd = func( p );
    //p += fd.xy*dt;                // gradient  move
    p += vec2( -fd.y, fd.x ) * dt;  // vorticity move     
}

vec3 bgtex( vec2 p ){
    vec2 sp = sin(p*10.0); 
    float w2 = 0.16;
    //float grid = w2/(w2 + dot(sp,sp));              // grid only nodes
    sp = w2/(w2 + sp*sp); float grid = dot(sp,sp);    // grid only lines
    //sp = w2/(w2 + sp*sp); float grid = sp.x*sp.y;   // grid with nodes, lines
    //return vec3( val );
    vec2 bg = sin(p*2.0)*0.5 + 0.5;
    //return vec3( bg.x, grid*0.5, bg.y ); 
    return vec3( bg.x, (bg.x+bg.y)*0.4, bg.y ) + vec3(grid*0.3); 
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
    //return vec4( vec3(fd.xy,0.0)*0.1 + vec3(fd.z*0.5+0.5), 0.0);
    //return vec4( vec3(fd.z*0.5+0.5), 0.0);
    return vec4( vec3(fd.xy,0.0)*0.1 + 0.5, 0.0);
}


void main( ){	
    
	vec2 p    = 2.0 * gl_FragCoord.xy / resolution.xy - 1.0;
	p *=6.0;	
	
	float val = func( p ).z;
	//float dt = 0.0625;
	float sint = sin(time*0.1);
	float wt2 = 0.4;
	float dt = (wt2/(wt2+sint*sint)-(wt2/(1.0+wt2)))*0.0625*32.0;
	vec4 clr = vec4(0.0,0.0,0.0,1.0);
	float t = 0.0;
	for(int i=0; i<32; i++){ 
	    move( p, dt ); 
	    //clr.xyz += bgtex( p )*vec3(t,0.5,1.0-t)*dt*4.0;
	    //clr.xyz += bgtex( p )*t*dt*4.0;
	    t += dt;
	};
    clr.xyz += bgtex( p ); 
    //clr.z = val;
    gl_FragColor = clr; 
    
    //vec4 clr = vec4( vec3(func(p).z+0.5), 1.0 );
    //vec4 clr = vec4( vec3(func(p).z+0.5), 1.0 );
    //gl_FragColor = checkFuncDeriv(p); 
}
