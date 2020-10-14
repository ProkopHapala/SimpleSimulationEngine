#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

void main(){
    vec4 c = textureLod( texture_1, fUV, 0 );
    //c.xy*=0.5; c.xy+=vec2(0.5);
	gl_FragColor= c;
	//gl_FragColor   = vec4( fUV, sin(fUV.x*10.0)*sin(fUV.y*10.0), 1.0 );
	//gl_FragColor = vec4(fUV.x,fUV.y,0.0,1.0);
	//gl_FragColor = vec4( fUV, 0.,1.);
}


