#version 330 core

smooth in vec4 gl_FragCoord;
in        vec2 gl_PointCoord;

out       vec4 gl_FragColor;

void main(){
    vec2  fUV    = 2.0*gl_PointCoord - vec2(1.0,1.0);
    float r2     = dot(fUV,fUV);
    
    //float alpha = (1-r2); alpha*=alpha;
    float alpha = 1/(1+r2*32.0);
    if( alpha > 0.05 ){
        gl_FragColor = vec4( 1-r2, 0.5, r2*8.0, alpha );
    }else{
        discard;
    }
}
