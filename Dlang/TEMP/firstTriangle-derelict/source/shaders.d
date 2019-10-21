module shaders;

string vertexShaderSrc = 
q{
    #version 330 core
    layout(location = 0) in vec3 inVertexPosition;
    //noperspective out vec3 world_pos;
    void main(){
        gl_Position = vec4(inVertexPosition, 1.0);
    }
};

string fragmentShaderSrc = 
q{
    #version 330 core
    smooth in vec4 gl_FragCoord;
    out vec3 color;
    void main(){
        //color = vec3(1,1,0);
        //color = sin( gl_FragCoord.xyz );
        color = vec3( sin( length(gl_FragCoord.xy) ), sin( gl_FragCoord.xy) );
        //color = vec3( sin( length(gl_FragCoord.xy) ) );
    }
};
