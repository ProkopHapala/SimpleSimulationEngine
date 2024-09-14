// The MIT License
// Copyright © 2014 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Signed distance to a triangle. Negative in the inside, positive in the outside.
// Note there's only one square root involved. The clamp(x,a,b) is really just
// max(a,min(b,x)). The sign(x) function is |x|/x. 


// List of other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and iquilezles.org/articles/distfunctions2d


// Other triangle functions:
//
// Distance:   https://www.shadertoy.com/view/XsXSz4
// Gradient:   https://www.shadertoy.com/view/tlVyWh
// Boundaries: https://www.shadertoy.com/view/tlKcDz


// signed distance to a 2D triangle
float sdTriangle( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 )
{
	vec2 e0 = p1 - p0;
	vec2 e1 = p2 - p1;
	vec2 e2 = p0 - p2;

	vec2 v0 = p - p0;
	vec2 v1 = p - p1;
	vec2 v2 = p - p2;

	vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );
	vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
	vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );
    
    float s = e0.x*e2.y - e0.y*e2.x;
    vec2 d = min( min( vec2( dot( pq0, pq0 ), s*(v0.x*e0.y-v0.y*e0.x) ),
                       vec2( dot( pq1, pq1 ), s*(v1.x*e1.y-v1.y*e1.x) )),
                       vec2( dot( pq2, pq2 ), s*(v2.x*e2.y-v2.y*e2.x) ));

	return -sqrt(d.x)*sign(d.y);
}

//---- everything below this line is NOT the triangle formula, and it's there just to produce the picture ---

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 p = (2.0*fragCoord.xy-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
	p *= 1.5;
    m *= 1.5;
    
    // animate
	vec2 v1 = vec2(1.4,1.0)*cos( iTime + vec2(0.0,2.00) + 0.0 );
	vec2 v2 = vec2(1.4,1.0)*cos( iTime + vec2(0.0,1.50) + 1.5 );
	vec2 v3 = vec2(1.4,1.0)*cos( iTime + vec2(0.0,3.00) + 4.0 );

    // distance
	float d = sdTriangle( p, v1, v2, v3 );
    
    // color
#if 1
    vec3 col = vec3(1.0) - sign(d)*vec3(0.1,0.4,0.7);
	col *= 1.0 - exp(-2.0*abs(d));
	col *= 0.8 + 0.2*cos(120.0*d);
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.02,abs(d)) );
#else
    vec3 col = (d<0.0) ? vec3(0.2,0.8,0.8) : vec3(0.8,0.2,0.8);
	col *= exp2(-2.0*abs(d));
    col = mix(col,vec3(1.20),exp2(-22.0*abs(d)));
	col = mix( col, vec3(0.0), 1.0-smoothstep(0.01,0.02,abs(d)) );
    col *= 1.0 + 0.5*smoothstep(-0.4,0.4,cos(100.0*d));
#endif

    // interactivity
    if( iMouse.z>0.001 )
    {
    d = sdTriangle(m, v1, v2, v3 );
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

	fragColor = vec4(col,1.0);
}