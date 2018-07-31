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

uniform  float txScale;


out vec4 gl_FragColor;

void main(){
    
    vec3  p     = fragUVW;
    
    
    /*
    float cover = 0.0;
    float t     = 0.0;
    vec3 dp     = fragUVWdir*0.01;
    for(int i=0; i<64; i++ ){
        float dens = textureLod( texture_1, p*txScale-vec3(0.5), 0 ).r;
        cover += dens;
        if(cover>=1.0) break;
        float dt = 1/(0.1+dens);
        //t += 0.5;
        t += dt;
        p += dp*dt; // TODO  non-uniform step need to modify integration scheme
    }
    //gl_FragColor = vec4(vec3(1-t*0.005),1.0);
    gl_FragColor = vec4(vec3(1.0,0.0,0.0),1-t*0.001);
    */
    
    /*
    vec3 dp     = fragUVWdir*0.1;
    float cover = 0.0;
    for(int i=0; i<16; i++ ){
        float dens = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        cover += dens;
        p += dp; // TODO  non-uniform step need to modify integration scheme
    }
    gl_FragColor = vec4(vec3(1.0,0.0,0.0),cover );
    */
    
    
    // self-shade - this kind of works but very costly
    /*
    vec3 dpp    = normalize(vec3(0.0,1.0,0.5))*0.1;
    vec3 dp     = fragUVWdir*0.1;
    float cover = 0.0;
    float light = 0.0;
    
    for(int i=0; i<16; i++ ){
        float dens = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        cover += dens;
        float enlight = 1.0;
        
        if(dens>0.05){ // march toward light
            vec3 pp  = p;
            for(int i=0; i<16; i++ ){
                pp += dpp;
                float dens2 = textureLod( texture_1, pp*txScale+vec3(0.5), 0 ).r;
                enlight *= (1-dens2*5.0);
            }
        }
        
        light += dens*enlight;
        if(cover>1.0) break;
        p += dp; // TODO  non-uniform step need to modify integration scheme
    }
    gl_FragColor = vec4(vec3(light),cover );
    */
    
    /*
    /// Self Shade cog - cast ray from density centre of mass 
    vec3  dp     = fragUVWdir*0.1;
    float cover  = 0.0;
    vec3  cog    = vec3(0.0);
    for(int i=0; i<64; i++ ){
        float dens = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        cover += dens;
        cog   += p*dens*max(0.0,1-cover);
        if(cover>1.0) break;
        p += dp; // TODO  non-uniform step need to modify integration scheme
    }

    dp    = normalize(vec3(0.0,1.0,0.5))*0.05;
    float light = 1.0;
    p = cog/cover;
    for(int i=0; i<32; i++ ){
        p += dp;
        float dens = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        light *= max(0.0,1-dens*0.5);
    }
    gl_FragColor = vec4(vec3(light*light),cover );
    */
    
    /*
    // efficient ray-marching
    // we need to modify texture as signed-distance-field
    //https://www.shadertoy.com/view/XsSBWc
    //float it = 0.0;
    float sum = 1.0;
    for(int i=0; i<64; i++){
        vec3 pos = ro + t*rd;
        float f = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        float dens = f*f;
        // --- nonuniform walk occluding
        float dt    = 0.05/(dens*8.0 + 1.0 );
        float da    = dens*dt*4.0;
        float wa    = 1.0-clrsum.a; wa=clamp(wa,0.0,1.0); 
        sum +=  da * wa;
        p += dp*dt;
        it+= (1.0/64.0);
    }
    gl_FragColor = vec4(vec3(sum),cover );
    */
    
    float trashold = 0.05;
    vec3 dp     = normalize(fragUVWdir);
    float sum = 0.0;
    float it = 0.0;
    for(int i=0; i<32; i++){
        float dist = textureLod( texture_1, p*txScale+vec3(0.5), 0 ).r;
        float dt = (0.0+max(0,dist-trashold) ); // overexted make sense in combination with regula-falsi
        it  +=1.0;
        sum += dt;
        p += dp*dt;
        if(dist<trashold) break;
        //if( sum>3.0 ) discard;
    }
    
    gl_FragColor = vec4(vec3(sum*0.3), 1.0 );  // ray length
    //gl_FragColor = vec4(vec3(it*0.025), 1.0 ); // iteration count

    vec4 rgba = textureLod( texture_1, p*txScale+vec3(0.5), 0 );
    //gl_FragColor = vec4( vec3(rgba.x),1.0);
}
