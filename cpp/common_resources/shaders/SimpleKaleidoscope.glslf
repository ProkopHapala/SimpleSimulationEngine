#version 330 core

in       vec2      fUV;
uniform  vec2      Const;    // julia constant

out vec4 gl_FragColor;

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = (fragCoord.xy-.5*iResolution.xy) * 7.2 / iResolution.y;

    float r = 1.0;
    float a = iTime*.1;
    float c = cos(a)*r;
    float s = sin(a)*r;
    for ( int i=0; i<32; i++ )
    {
    	uv = abs(uv);
        uv -= .25;
        uv = uv*c + s*uv.yx*vec2(1,-1);
    }
        
    fragColor = .5+.5*sin(iTime+vec4(13,17,23,1)*texture( iChannel0, uv*vec2(1,-1)+.5, -1.0 ));
}
