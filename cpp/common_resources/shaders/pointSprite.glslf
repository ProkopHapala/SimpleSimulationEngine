#version 330 core

// gl_PointCoord is a built-in input variable for fragment shaders when rendering points.
// No explicit 'in' declaration is needed in core profile.

out vec4 fragColor; // Define a custom output color variable

void main(){
    vec2  fUV    = 2.0*gl_PointCoord - vec2(1.0,1.0);
    float r2     = dot(fUV,fUV);
    
    float alpha = 1/(1+r2*32.0);
    if( alpha > 0.05 ){
        fragColor = vec4( 1-r2, 0.5, r2*8.0, alpha ); // Assign to the custom output variable
    }else{
        discard;
    }
}
