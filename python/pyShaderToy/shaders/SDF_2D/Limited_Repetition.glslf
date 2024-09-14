// The MIT License
// Copyright © 2019 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// This shader shows the correct way to compute a number of
// copies of an object (right), and the usual but incorrect
// way to do it (left). The incorrect way is to do infinite
// repetition and then clip it with a box.
//
// More info: https://iquilezles.org/articles/sdfrepetition/


// Create multiple copies of an object - https://iquilezles.org/articles/sdfrepetition/
vec2 opRepLim( in vec2 p, in float s, in vec2 lima, in vec2 limb )
{
    return p-s*clamp(round(p/s),lima,limb);
}

// Create infinite copies of an object - https://iquilezles.org/articles/sdfrepetition/
vec2 opRep( in vec2 p, in float s )
{
    return p-s*round(p/s);
    //return mod(p+s*0.5,s)-s*0.5;
}

// https://iquilezles.org/articles/distfunctions
float sdBox( in vec2 p, in vec2 b ) 
{
    vec2 q = abs(p) - b;
    return min(max(q.x,q.y),0.0) + length(max(q,0.0));
}

//-----------------------------

float map( in vec2 p )
{
    float d;
    if( p.x>0.0 ) // correct way to do it
    {
    	vec2 q = p*6.0 - vec2(5.0,0.0);
        vec2 r = opRepLim(q,2.0,vec2(-1,-2),vec2(1,2));
        d = sdBox( r, vec2(0.4,0.2) ) -  0.1;
    }
    else         // incorrect way to do it
    {
    	vec2 q = p*6.0 + vec2(5.0,0.0);
        vec2 r = opRep(q,2.0);
        d = sdBox( r, vec2(0.4,0.2) ) -  0.1;
        d = max( d, sdBox( q, 2.0*vec2(1,2)+0.5 ) );
    }
    return d/6.0;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
	
    // sdf
    float d = map(p);
    
    // colorize
    vec3 col = (d>0.0) ? vec3(0.9,0.6,0.3) : vec3(0.65,0.85,1.0);
	col *= 1.0 - exp(-24.0*abs(d));
	col *= 0.8 + 0.2*cos(240.0*d);
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.01,abs(d)) );

    if( iMouse.z>0.001 )
    {
    d = map( m );
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, max(abs(length(p-m)-abs(d))-0.0025,-sign(m.x)*p.x)));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015) );
    }

    col *= smoothstep(0.005,0.010,abs(p.x));
    
	fragColor = vec4(col,1.0);
}