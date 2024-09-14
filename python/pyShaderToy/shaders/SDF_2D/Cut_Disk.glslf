// The MIT License
// Copyright © 2022 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// Signed distance to a disk that's been clipped by a line
//
// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and iquilezles.org/articles/distfunctions2d



// r=radius, h=height
float sdCutDisk( in vec2 p, in float r, in float h )
{
    float w = sqrt(r*r-h*h); // constant for a given shape
    
    p.x = abs(p.x);
    
    // select circle or segment
    float s = max( (h-r)*p.x*p.x+w*w*(h+r-2.0*p.y), h*p.x-w*p.y );

    return (s<0.0) ? length(p)-r :        // circle
           (p.x<w) ? h - p.y     :        // segment line
                     length(p-vec2(w,h)); // segment corner
}



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // normalized pixel coordinates
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
   
    // animation
    float ra = 0.75;
    float he = ra*clamp(cos(iTime*0.8),-0.999999,0.999999);
   
    // distance
    float d = sdCutDisk(p,ra,he);
   
    // coloring
    vec3 col = (d>0.0) ? vec3(0.9,0.6,0.3) : vec3(0.5,0.85,1.0);
    col *= 1.0 - exp(-7.0*abs(d));
    col *= 0.8 + 0.2*cos(128.0*abs(d));
    col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.015,abs(d)) );

    // interactivity
    if( iMouse.z>0.001 )
    {
    d = sdCutDisk(m,ra,he);
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

	fragColor = vec4(col, 1.0);
}