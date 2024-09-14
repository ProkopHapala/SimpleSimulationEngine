// The MIT License
// Copyright © 2019 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// Distance to a U shape, based on https://www.shadertoy.com/view/tsjGzD
//
// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and https://iquilezles.org/articles/distfunctions2d


float sdHorseshoe( in vec2 p, in vec2 c, in float r, in float le, float th )
{
    p.x = abs(p.x);
    float l = length(p);
    p = mat2(-c.x, c.y, 
              c.y, c.x)*p;
    p = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
             (p.x>0.0)?p.y:l );
    p = vec2(p.x-le,abs(p.y-r)-th);
    return length(max(p,0.0)) + min(0.0,max(p.x,p.y));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // normalized pixel coordinates
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
    
    // animation
    float t =            3.14* (0.5+0.5*cos(iTime*0.5));
    vec2  w = vec2(0.750,0.25)*(0.5+0.5*cos(iTime*vec2(0.7,1.1)+vec2(0.0,3.0)));
    
    // distance
    float d = sdHorseshoe(p-vec2(0.0,-0.1),vec2(cos(t),sin(t)), 0.5, w.x, w.y);
        
    // coloring
    vec3 col = vec3(1.0) - sign(d)*vec3(0.1,0.4,0.7);
	col *= 1.0 - exp(-6.0*abs(d));
	col *= 0.8 + 0.2*cos(120.0*abs(d));
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.02,abs(d)) );

    if( iMouse.z>0.001 )
    {
    d = sdHorseshoe(m-vec2(0.0,-0.1),vec2(cos(t),sin(t)), 0.5, w.x, w.y);
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

	fragColor = vec4(col, 1.0);
}