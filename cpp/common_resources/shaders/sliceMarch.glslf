#version 330 core


//smooth in vec3 fragUVW;
//smooth in vec3 fragUVWstep;
noperspective in vec3 fragUVW;
noperspective in vec3 fragUVWdir;
uniform sampler3D texture_1; 
uniform sampler3D texture_noise;
uniform  float txScale;

out vec4 gl_FragColor;


void main(){

    vec3 gridStep    = vec3( 1.0/4.0, 1.0/4.0, 1.0/4.0 );
    vec3 invGridStep = vec3( 4.0, 4.0, 4.0 );

    vec3 p = fragUVW*txScale + vec3(0.5);
    vec3 dp = normalize(fragUVWdir);

    vec3 invDp = 1.0/dp;

    //vec3 dg = dp*invGridStep;

    vec4  color = vec4(0.0);

    for(int i=0; i<16; i++){
       vec3 ts  = p*invDp;
       vec3 tds = fract( ts ); // how much remains to next slice in each x,y,z direction ?
       vec3 tis = floor( ts ); // slice index in x,y,z direction
       //float t = min( min(tds.x,tds.y), tds.z );
       vec2 uv;
       if( tds.x<tds.y ){
            if(tds.x<tds.z){ uv=p.zy + vec2(0,tids.x); p += dp*tds.x; } // x=min
            else           { uv=p.xy + vec2(2,tids.z); p += dp*tds.z; } // z=min
       }else{
            if(tds.y<tds.z){ uv=p.xz + vec2(1,tids.y); p += dp*tds.y; } // y=min
            else           { uv=p.xy + vec2(2,tids.z); p += dp*tds.z; } // z=min
       }
       float val = texture( texture_1, uv, );
       color += val;
       if(color.a>1.0) break;
    }
    gl_FragColor = color;
}




