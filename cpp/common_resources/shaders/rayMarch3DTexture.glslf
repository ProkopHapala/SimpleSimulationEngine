#version 330 core

//smooth in vec3 fragUVW;
//smooth in vec3 fragUVWstep;
noperspective in vec3 fragUVW;
noperspective in vec3 fragUVWdir;
uniform sampler3D texture_1; 

out vec4 gl_FragColor;

void main(){
    vec3  color = vec3(0.0);
    float cover = 0.0;
    vec3  p     = fragUVW;
    vec3 dp     = fragUVWdir*0.01;
    for(int i=0; i<64; i++ ){
        vec4 tex = textureLod( texture_1, p, 0 );
        color += tex.rgb*tex.a;
        cover += tex.a;
        if(cover>=1.0) break;
        p += dp/(0.1+tex.a); // TODO  non-uniform step need to modify integration scheme
        //p+=dp;
    }
    gl_FragColor = vec4(color,cover);
    //gl_FragColor = vec4(color,1.0);
    //gl_FragColor = vec4(fragUVWdir+vec3(0.5,0.5,0.5),1.0);
}
