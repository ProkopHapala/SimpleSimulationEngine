#version 330 core


///  see this:
//   https://www.shadertoy.com/view/Md2BWc
// Approximation of the real function, computed with Maple
#define DENS 3.          // tau.rho at the center
#define rad 0.9     	 // sphere radius
vec3  sphericalTransmittanceGradient(vec2 L, float r, float h, float z){
	float Lx=L.x;
	float Ly=L.y;
	float Y = (r -h);
	float xmax = sqrt(2.0*r*h - h*h);
	float b = xmax;
	float a1 = (-xmax*0.7);
	if (DENS < 2.)
		a1 = 0.0;
	float a12 = a1*a1;float a13 = a12*a1;float a14 = a12*a12; 
	float Lx2 = Lx*Lx;float Lx3 = Lx2*Lx;float Lx4 = Lx3*Lx; float Lx5 = Lx4*Lx;float Lx6 = Lx5*Lx;
	float Ly2 = Ly*Ly;float Ly3 = Ly2*Ly;float Ly4 = Ly2*Ly2;float Ly5 = Ly4*Ly;float Ly6 = Ly5*Ly;
	float xmax3 = xmax*xmax*xmax;
	float Y2 = Y*Y;float Y3 = Y2*Y;float Y4 = Y2*Y2;
	float r2 = r*r;float r4 = r2*r2;
	float R2 = rad*rad;
	float H2 = z*z;
	float S = sqrt(a12*Lx2+Y2*Ly2-a12+2.*a1*Lx*Y*Ly-Y2+r2);
	float c1 = S*xmax3+5.*Lx2*r2*Y2*Ly2-3.*R2*Y2*Ly2+3.*H2*Y2*Ly2-5.*Lx2*Y4*Ly2
		-2.*Lx2*Y2*r2+5.*Ly4*Y2*r2-8.*Y2*Ly2*r2+4.*Lx2*Y4*Ly4-2.*S*a13
		-21.*S*a12*Lx2*Y*Ly+12.*S*Ly3*a12*Lx2*Y+12.*S*Lx4*a12*Y*Ly
		-3.*S*Lx2*Y*Ly*r2-2.*Ly2*a14+22.*Lx4*a14-8.*Lx6*a14-20.*a14*Lx2
		-3.*a12*r2+3.*Y2*a12-4.*Y2*Ly2*a12+Ly4*Y2*a12-8.*Ly2*a14*Lx4
		+4.*Lx4*a12*Y2-7.*Y2*a12*Lx2+10.*Ly2*a14*Lx2+Ly2*a12*r2-4.*Lx4*a12*r2
		+7.*a12*Lx2*r2+6.*a14-20.*Ly3*a13*Lx3*Y-12.*Ly4*a12*Lx2*Y2+11.*Ly3*a13*Lx*Y
		-20.*Lx5*a13*Y*Ly-12.*Lx4*a12*Y2*Ly2+41.*Lx3*a13*Y*Ly+23.*Lx2*Y2*Ly2*a12
		-21.*a13*Lx*Y*Ly+4.*a1*Lx3*Y3*Ly3-7.*a1*Y3*Ly3*Lx+3.*a1*Y3*Lx*Ly
		+4.*a1*Ly5*Y3*Lx-a1*Lx3*Y3*Ly-4.*Ly2*a12*Lx2*r2-6.*S*Y3*Ly+9.*S*Y3*Ly3
		+3.*S*H2*xmax+3.*S*Y2*xmax+3.*R2*Y2-3.*R2*r2-3.*H2*Y2+3.*H2*r2+10.*Y4*Ly2
		+3.*Y2*r2+Lx2*Y4+4.*Ly6*Y4-11.*Ly4*Y4+Ly2*r4+Lx2*r4-3.*Y4-4.*S*Ly5*Y3
		+8.*S*Lx5*a13-3.*S*R2*xmax-18.*S*a13*Lx3+12.*S*a13*Lx+3.*S*R2*Y*Ly
		-6.*S*Ly2*a13*Lx+8.*S*Ly2*a13*Lx3+6.*S*a12*Y*Ly-3.*S*H2*Y*Ly+3.*S*Lx2*Y3*Ly
		+3.*S*Y*Ly*r2-4.*S*Lx2*Y3*Ly3-3.*S*Ly3*Y*a12-3.*S*Ly3*Y*r2-3.*a1*R2*Lx*Y*Ly
		+3.*a1*H2*Lx*Y*Ly+a1*Ly3*Y*Lx*r2+a1*Lx3*Y*Ly*r2;	
	c1 *= (1./3.)*DENS/(S*R2);	
	float c2 = Y2*S-4.*Ly4*Y2*Lx*S+2.*Ly3*Y*S*a1-4.*Ly2*a12*Lx3*S
		+3.*Ly2*a12*Lx*S-8.*Lx4*a1*Y*Ly*S+14.*Lx2*Y*Ly*S*a1-3.*a13
		-4.*Y*Ly*S*a1-Ly2*S*Lx*r2-4.*Lx3*Y2*Ly2*S+7.*Y2*Ly2*Lx*S
		+9.*Lx3*a12*S+R2*Lx*S-2.*Y2*Lx*S-Lx3*S*r2+Lx*S*r2-H2*Lx*S
		-6.*a12*Lx*S-4.*Lx5*a12*S+Lx3*S*Y2-R2*S+a12*S-8.*Ly3*a1*Lx2*Y*S
		+12.*Ly3*a12*Lx3*Y+12.*Ly4*a1*Lx2*Y2-7.*Ly3*a12*Lx*Y+12.*Lx5*a12*Y*Ly
		+12.*Lx4*a1*Y2*Ly2-25.*Lx3*a12*Y*Ly-23.*Lx2*Y2*Ly2*a1+13.*a12*Lx*Y*Ly
		-R2*Lx*Y*Ly+H2*Lx*Y*Ly+5.*Y2*Ly2*a1-2.*Ly4*Y2*a1+4.*Ly2*a13*Lx4
		-3.*Lx4*a1*Y2+4.*Lx3*Y3*Ly3+6.*Y2*a1*Lx2-9.*Y3*Ly3*Lx
		+5.*Y3*Lx*Ly-R2*a1*Lx2+H2*a1*Lx2-5.*Ly2*a13*Lx2+4.*Ly5*Y3*Lx-3.*Lx3*Y3*Ly
		-Ly2*a1*r2+3.*Lx4*a1*r2-5.*a1*Lx2*r2+2.*a1*r2-11.*Lx4*a13+4.*Lx6*a13
		+10.*a13*Lx2+H2*S+R2*a1-3.*Y2*a1-H2*a1+Ly2*a13+3.*Ly2*a1*Lx2*r2
		+3.*Ly3*Y*Lx*r2+3.*Lx3*Y*Ly*r2-4.*Lx*Y*Ly*r2;
	c2 *= DENS/(R2*S);
	if (abs(c2) < 0.1)
		c2 = 0.1; // arbitraire
	float EX1 = exp(c1-c2*xmax);
	float EX2 = exp(c1+c2*xmax);
	float res = -2.*EX1+EX1*c2*c2*R2-EX1*c2*c2*Y2-EX1*c2*c2*H2
		-2.*EX1*c2*xmax-EX1*xmax*xmax*c2*c2+2.*EX2-EX2*c2*c2*R2+EX2*c2*c2*Y2+EX2*c2*c2*H2
		-2.*EX2*c2*xmax+EX2*xmax*xmax*c2*c2;
	res *= -DENS/(rad*rad*c2*c2*c2);
	return vec3(res);
}


//smooth in vec3 fragUVW;
//smooth in vec3 fragUVWstep;
noperspective in vec3 fragUVW;
noperspective in vec3 fragUVWdir;
uniform sampler3D texture_1; 
uniform sampler3D texture_noise;

uniform  float txScale;

out vec4 gl_FragColor;


const float safeFator = 1.0;

/*
float rayMarchDistTexture( inout vec3 p, vec3 dp, float iso, float overshoot ){
    //float t = 0.0;
    float dist;
    for(int i=0; i<32; i++){
        dist = textureLod( texture_1, p+vec3(0.5), 0 ).r;
        dist -= trashold;
        float dt = overshoot+safeFato*max(0.0,dist); // overexted make sense in combination with regula-falsi
        if(dist<0.0) break;
        //t += dt;
        p += dp*dt;
    }
    return dist;
}
*/

vec3 bump( vec3 p, float freq ){
 return vec3(
        textureLod( texture_noise, p*freq+vec3(0.5,0.0,0.0), 0 ).r-0.5,
        textureLod( texture_noise, p*freq+vec3(0.0,0.5,0.0), 0 ).r-0.5,
        textureLod( texture_noise, p*freq+vec3(0.0,0.0,0.5), 0 ).r-0.5
    );
}

float distField( vec3 p, float iso){
    return textureLod( texture_1, p, 0 ).r - iso;
    //return textureLod( texture_1, p + bump( p, 3.0 )*0.04 , 0 ).r - iso;
    //return textureLod( texture_1, p, 0 ).r - iso*( 1 + 0.5*textureLod( texture_noise, p*8.0, 0 ).r );
    //float fbig  = textureLod( texture_1, p, 0 ).r - iso;
    //float ffine = textureLod( texture_noise, p*5.0, 0 ).r;
    //return (fbig>0)? fbig : ffine*fbig;
}

float distFieldFine( vec3 p, float iso){
    return textureLod( texture_1, p, 0 ).r - iso;
    //return textureLod( texture_1, p + bump( p, 3.0 )*0.02 , 0 ).r;
    //return textureLod( texture_1, p, 0 ).r - iso*( 1 + 0.5*textureLod( texture_noise, p*8.0, 0 ).r );
    //float fbig  = textureLod( texture_1, p, 0 ).r - iso;
    //float ffine = textureLod( texture_noise, p*5.0, 0 ).r;
    //return (fbig>0)? fbig : ffine*fbig;
}

vec3 getGradient( vec3 p, float d, float iso){
    // TODO - perhaps we can use hardware derivatives?
    //float y = distField(p,iso); // should be zero? why to compute it again ?
    //return vec3( distField( p+vec3(d,0.0,0.0), iso) - y,
    //             distField( p+vec3(0.0,d,0.0), iso) - y,
    //             distField( p+vec3(0.0,0.0,d), iso) - y  ) / d;
    float y = distFieldFine(p,iso); // should be zero? why to compute it again ?
    return vec3( distFieldFine( p+vec3(d,0.0,0.0), iso) - y,
                 distFieldFine( p+vec3(0.0,d,0.0), iso) - y,
                 distFieldFine( p+vec3(0.0,0.0,d), iso) - y  ) / d;
}

vec3 getGradient2( vec3 p, float d, float iso){
    return vec3( distFieldFine( p+vec3(d,0.0,0.0), iso) - distFieldFine( p+vec3(-d,0.0,0.0), iso),
                 distFieldFine( p+vec3(0.0,d,0.0), iso) - distFieldFine( p+vec3(0.0,-d,0.0), iso),
                 distFieldFine( p+vec3(0.0,0.0,d), iso) - distFieldFine( p+vec3(0.0,0.0,-d), iso) ) / (d*2);
}



float rayMarchDistTexture2( inout vec3 p, vec3 dp, float iso, float overshoot ){
    float t = 0.0;
    float y,oy;
    float dt;
    for(int i=0; i<32; i++){
        y  = distField(p, iso);
        dt = overshoot+safeFator*max(0.0,y); // overexted make sense in combination with regula-falsi
        if(y<0.0) break;
        if(t>1.7) return t;
        t += dt;
        oy = y;
        p += dp*dt;
    }
    // linear regresion ( 1 iter of regula-falsa .. we can do more )
    float f = -oy/(y-oy);
    dt *= (f-1.0);
    t += dt;
    p += dp*dt;
    return t;
}

float rayMarchDistTexture2Leafs( inout vec3 p, vec3 dp, float iso, float overshoot ){
    float t = 0.0;
    float y,oy;
    float dt;
    for(int i=0; i<16; i++){
        y  = distField(p,iso);
        dt = overshoot+safeFator*max(0.0,y); // overexted make sense in combination with regula-falsi
        if( y<0 ){
            float mask = textureLod( texture_noise, p*5.0, 0 ).r;
            if (mask>0.5) return t;
        }
        t += dt;
        oy = y;
        p += dp*dt;
    }
    return t;
}

float rayMarchDistTexture2Bumpy( inout vec3 p, vec3 dp, float iso, float overshoot ){
    float t = 0.0;
    float y,oy;
    float dt;
    for(int i=0; i<16; i++){
        y  = distField(p, iso);
        if( y<0 ){
            float mask = 0.5 - textureLod( texture_noise, p*5.0, 0 ).r;
            //y += 0.5 - maks;
            //if( y<0.0 ) break;
            dt = 0.01;
            if(mask<0) break;
        }else{
            dt = overshoot+safeFator*max(0.0,y);
        }
        t += dt;
        oy = y;
        p += dp*dt;
    }
    // linear regresion ( 1 iter of regula-falsa .. we can do more )
    float f = -oy/(y-oy);
    dt *= (f-1.0);
    t += dt;
    p += dp*dt;
    return t;
}

float rayMarchShadow( inout vec3 p, vec3 dp, float iso, float dt0, float dtfac, float dens2 ){
    float t = 0.0;
    float y;
    float dt    = dt0;
    float light = 1.0;
    for(int i=0; i<16; i++){
        y      = distField(p, iso);
        light *= clamp(1 + y*dens2, 0.0, 1.0 );
        dt *= dtfac;
        t  += dt;
        p  += dp*dt;
        if( (p.x<0.1)||(p.x>0.9)||(p.y<0.1)||(p.y>0.9)||(p.z<0.1)||(p.z>0.9) ) return light;
        //if(t>1.7) return light;
    }
    return light;
}

void main(){

    vec3 lightDir = normalize(vec3(0.0,1.0,-0.5));

    vec3 p = fragUVW*txScale + vec3(0.5);
    vec3 dp = normalize(fragUVWdir);
    p+=dp*0.1;

    float iso = 0.06;

    float t = rayMarchDistTexture2( p, dp, 0.04, 0.1/16.0 );
    //float t = rayMarchDistTexture2Leafs( p, dp, 0.06, 0.5/16.0 );
    //float t = rayMarchDistTexture2Bumpy( p, dp, 0.06, 0.5/16.0 );

    /*
    float mask;
    mask = textureLod( texture_noise, p*5.0, 0 ).r;
    if(mask<0.5) float t = rayMarchDistTexture2( p, dp, 0.05, 0.2/16.0 );
    mask = textureLod( texture_noise, p*5.0, 0 ).r;
    if(mask<0.5) float t = rayMarchDistTexture2( p, dp, 0.04, 0.2/16.0 );
    mask = textureLod( texture_noise, p*5.0, 0 ).r;
    if(mask<0.5) float t = rayMarchDistTexture2( p, dp, 0.02, 0.2/16.0 );
    */

    if( (p.x<0.1)||(p.x>0.9)||(p.y<0.1)||(p.y>0.9)||(p.z<0.1)||(p.z>0.9) ) discard;

    //p += bump(p, 5.0) * (0.5/16.0);
    //p += dp * normalize( grad );

    //vec3  grad   = getGradient( p, 0.01, iso );
    vec3  grad   = getGradient( p, 0.5/16.0, iso );     // larger derivative step "d" leads to smooth shading, smaller makes flat shading
    //vec3  grad   = getGradient2( p, 0.25/16.0, iso );

    //float d = (0.25/16.0);
    //vec3  grad = ( getGradient( p+vec3(d,0.0,0.0), 0.01, iso ) +
    //               getGradient( p+vec3(0.0,d,0.0), 0.01, iso ) +
    //               getGradient( p+vec3(0.0,0.0,d), 0.01, iso ) ) * 0.3;

    //grad        += bump(p, 5.0);
    vec3  normal = normalize( grad );

    float light  = 1.0;
    light = max(0.0,dot( lightDir, normal ))*0.8 + 0.2;
    //light = textureLod( texture_noise, p*5.0, 0 ).r;
    //light     =  rayMarchShadow( p, vec3(0.0,1.0,0.0), 0.1, 0.05, 1.0, 2.5 );
    gl_FragColor = vec4( (normal*0.5 + 0.5)*light, 1.0 );
    //gl_FragColor = vec4( vec3(light), 1.0 );

    //gl_FragColor = vec4(vec3(t*0.3), 1.0 );  // ray length
    //gl_FragColor = vec4(vec3(it*0.025), 1.0 ); // iteration count

    //vec4 rgba = textureLod( texture_1, p*txScale+vec3(0.5), 0 );
    //gl_FragColor = vec4( vec3(rgba.x),1.0);
}




