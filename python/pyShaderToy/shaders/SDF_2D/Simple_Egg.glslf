// The MIT License
// Copyright © 2020 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// A simple generaliztion of sylvain69780's shader "Moss's EGG":
// https://www.shadertoy.com/view/wsBBR3
//
// This is also a special case of the Vesica primitive, in particular
// this egg is just half a vesica with a circle attached at its end
// https://www.shadertoy.com/view/XtVfRW

// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and https://iquilezles.org/articles/distfunctions2d


float sdEgg( in vec2 p, in float he, in float ra, in float rb )
{
    float ce = 0.5*(he*he-(ra-rb)*(ra-rb))/(ra-rb);

    p.x = abs(p.x);

    if( p.y<0.0 )             return length(p)-ra;
    if( p.y*ce-p.x*he>he*ce ) return length(vec2(p.x,p.y-he))-rb;
                              return length(vec2(p.x+ce,p.y))-(ce+ra);
}

float sdLineV( in vec2 p, in float b ) { p.y -= clamp( p.y, 0.0, b ); return length( p ); }

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // normalized pixel coordinates
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
    float px = 2.0/iResolution.y;
    
    p.y += 0.35;
    m.y += 0.35;

    // animation
    vec2 cen = vec2(0.0,0.0);
    float he = 0.7 + 0.3*cos(iTime*1.0+0.0);
    float ra = 0.5;
    float rb = 0.2;
    float al = smoothstep( -0.5, 0.5,sin(iTime+0.1) );
        
    // distance
    float d = sdEgg(p-cen, he, ra, rb);
    
    // coloring
    vec3 col = (d>0.0) ? vec3(0.9,0.6,0.3) : vec3(0.65,0.85,1.0);
	col *= 1.0 - exp(-7.0*abs(d));
	col *= 0.8 + 0.2*cos(90.0*abs(d));
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,1.5*px,abs(d)-0.005) );

    // draw primitive parameters
    d = abs( sdLineV(p-cen,he));
    d = min( d, abs(length(p)-ra) );
    d = min( d, abs(length(p-vec2(0.0,he))-rb) );
    d -= 0.003;
	col = mix( col, vec3(0.0,1.0,1.0), al*(1.0-smoothstep(0.0,1.5*px,d)) );

    if( iMouse.z>0.001 )
    {
    d = sdEgg(m-cen, he, ra, rb);
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

	fragColor = vec4(col, 1.0);
}