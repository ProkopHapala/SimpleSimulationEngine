// The MIT License
// Copyright © 2022 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// Signed distance to a staircase
//
// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and iquilezles.org/articles/distfunctions2d


float dot2( in vec2 v ) { return dot(v,v); }

float sdSquareStairs( in vec2 p, in float s, in float n )
{
    // constant for a given shape
    const float kS2 = sqrt(2.0);
    float w = 2.0*n+1.0;
    
    // pixel dependent computations
    p = vec2( abs(p.y+p.x), p.y-p.x ) * (0.5/s);

    float x1 = p.x-w;
    float x2 = abs(p.x-2.0*min(round(p.x/2.0),n))-1.0;
    
    float d1 = dot2( vec2(x1, p.y) + clamp(0.5*(-x1-p.y), 0.0, w  ) );
    float d2 = dot2( vec2(x2,-p.y) + clamp(0.5*(-x2+p.y), 0.0, 1.0) );

    return sqrt(min(d1,d2)) *
           sign(max(x1-p.y,(x2+p.y)*kS2)) *
           s*kS2;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // normalized pixel coordinates
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
 
    // animate
    float w = 1.0/8.0;
    float n = floor( 3.95*(0.5 + 0.5*cos(iTime*3.0)) );
 
    // distance
    float d = sdSquareStairs(p,w,n);
   
    // coloring
    vec3 col = (d>0.0) ? vec3(0.9,0.6,0.3) : vec3(0.65,0.85,1.0);
    col *= 1.0 - exp(-7.0*abs(d));
    col *= 0.8 + 0.2*cos(160.0*abs(d));
    col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.015,abs(d)) );

    // interactivity
    if( iMouse.z>0.001 )
    {
    d = sdSquareStairs(m,w,n);
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

	fragColor = vec4(col, 1.0);
}