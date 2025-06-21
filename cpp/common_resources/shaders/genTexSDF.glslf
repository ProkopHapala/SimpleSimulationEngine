#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

//out vec4 gl_FragColor;
out vec4 fragColor;

//uniform ivec2 n;
//#define msub 8
uniform int    msub;
uniform float  lpix;

void main(){
    float h =lpix*msub*0.5;
    float ih=1.0/h;
    float rOn  = 1.0;
    float rOff = 1.0;
    for(int ix=0; ix<msub; ix++ ){
        for(int iy=0; iy<msub; iy++ ){
            vec2   d  = vec2( ix, iy )*lpix-h;
            float  r  = length(d)*ih*0.70710678118;
            //float  r  = length(d)*ih*0.25;
            //float  r  = length(d)*ih*2.00;
            if(r<1.0){
                vec4 mask = texture( texture_1, fUV+d  );
                if(mask.r>0.5){ rOff = min(r,rOff); }
                else          { rOn  = min(r,rOn ); }
            }
            //rOn += texture( texture_1, fUV+d  ).r;
        }
    }
    float val = (rOn-rOff)*0.5 + 0.5;
    //float val = rOn/(msub*msub);
    //gl_FragColor = vec4( val, 0.0,0.0,1.0 );
    fragColor = vec4( val );
}


