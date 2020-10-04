#version 330 core

in       vec2      fUV;
uniform  vec2      Const;    // julia constant

out vec4 gl_FragColor;

vec2 mul_cmplx(vec2 a, vec2 b){
    return vec2( 
        a.x*b.x - a.y*b.y,
        a.x*b.y + a.y*b.x
    );
}

float juliaSet( vec2 z, vec2 c ){
    float ret = 0;
    for(int i=1;i<=256;i++){
        //z = mul_cmplx(z,z) + c;                            // Julia z^2 (standard)
        //z = mul_cmplx( mul_cmplx(z,z), z ) + c;              // Julia z^3
        //z = mul_cmplx( mul_cmplx(z,z), mul_cmplx(z,z) ) + c; // Julia z^4
        z = mul_cmplx( mul_cmplx( mul_cmplx(z,z), mul_cmplx(z,z) ), z ) + c; // Julia z^5
        if(dot(z,z)>4.0) break;
        ret+=1.0;
    }
    return ret;
}
 

void main(){
    //gl_FragColor   = textureLod( texture_1, fUV, 0 );
    vec2 z = (fUV*2.)+vec2(-2.0,-2.0);
    //vec2 Const = vec2(-0.5,0.3);
    float c = log(juliaSet( z, Const ) ) *0.2;
    gl_FragColor = vec4(c,c,c,1.);
    //gl_FragColor = vec4();
    //gl_FragColor = vec4(1.,1.,1.,1.);
}


