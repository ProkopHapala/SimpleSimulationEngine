#version 330 core

noperspective in vec3 fragUVW;
noperspective in vec3 fragUVWdir;
uniform sampler2D texture_1; 
//uniform sampler2D texture_noise;
uniform  float txScale;

out vec4 gl_FragColor;


void main(){

    vec3 gridStep    = vec3( 1.0/4.0, 1.0/4.0, 1.0/4.0 );
    vec3 invGridStep = vec3( 4.0, 4.0, 4.0 );

    vec3 p = fragUVW*txScale + vec3(0.5);
    vec3 dp = normalize(fragUVWdir)*0.25;
    vec3 invDp = 1.0/dp;

    vec2 duv    = vec2(4.0,3.0);
    vec2 invDuv = 1.0/duv;
    
    //vec3 dg = dp*invGridStep;

    
    vec4  color = vec4(0.0);
    for(int i=0; i<4; i++){
       vec3 ts  = p*invDp;
       vec3 tds = 1-fract( ts ); // how much remains to next slice in each x,y,z direction ?
       vec3 tis = floor( ts ); // slice index in x,y,z direction
       //float t = min( min(tds.x,tds.y), tds.z );
       vec2 uv;
       //uv=clamp(p.xy,0.01,0.99)*invDuv + vec2(tis.z*invDuv.x,invDuv.y*2);  p += dp*(tds.z+0.0001);
       uv=clamp(p.zy,0.01,0.99)*invDuv + vec2(tis.x*invDuv.x,0         ); p += dp*(tds.x+0.0001);
       /*
       if( tds.x<tds.y ){
            if(tds.x<tds.z){ uv=clamp(p.zy,0.01,0.99)*invDuv + vec2(tis.x*invDuv.x,0         ); p += dp*(tds.x+0.0001); } // x=min
            else           { uv=clamp(p.xy,0.01,0.99)*invDuv + vec2(tis.z*invDuv.x,invDuv.y*2); p += dp*(tds.z+0.0001); } // z=min
       }else{
            if(tds.y<tds.z){ uv=clamp(p.xz,0.01,0.99)*invDuv + vec2(tis.y*invDuv.x,invDuv.y  ); p += dp*(tds.y+0.0001); } // y=min
            else           { uv=clamp(p.xy,0.01,0.99)*invDuv + vec2(tis.z*invDuv.x,invDuv.y*2); p += dp*(tds.z+0.0001); } // z=min
       }
       */
       color += textureLod( texture_1, uv, 0 );
    }
    color /= 4.0;
    
    
    /*
    vec4 color = vec4(0.0);
    for(int i=0; i<4; i++){
       vec3 ts  = p*invDp;
       vec3 tds = fract( ts );
       vec3 tis = floor( ts );
       vec2 uv  = clamp(p.xy,0.01,0.99)*invDuv + vec2(tis.z*invDuv.x,0);
       p       += dp*( (1-tds.z)+0.0001);
       color += textureLod( texture_1, uv, 0 );
    }
    color /= 4.0;
    */
    
    //gl_FragColor = vec4( vec3(0.0), cover );
    //vec4 color = textureLod( texture_1, p.xy, 0 );
    
    gl_FragColor = color;
}




