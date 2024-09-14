// https://www.shadertoy.com/view/lfX3RX adaptive stairs, 2023 jt
// based on https://www.shadertoy.com/view/lclGz2 thales circle stairs
// based on https://www.shadertoy.com/view/XcX3zX thales circle box

// Stairs with adaptive steps.
// Use mouse to change slope.

// TODO: How much trigonometry can be removed?
// TODO: Implement thickness for triangular variant.

// tags: sdf, circle, distance, box, angle, euclidean, stairs, signed, adaptive, exact, inscribed, thales

// The MIT License
// Copyright (c) 2023 Jakob Thomsen
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#define pi 3.1415926

float box(vec2 p, vec2 b) // https://iquilezles.org/articles/distfunctions2d/
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
}

float hseg(vec2 p)
{
    return length(p - vec2(clamp(p.x,-1.0,+1.0),0));
}

float adaptive_stairs(vec2 p, float n) // https://www.shadertoy.com/view/lclGz2 (jt)
{
    float mu = pi/6.0*iMouse.x/iResolution.x;
    if(all(lessThan(iMouse.xy,vec2(10.0)))) mu = (0.5+0.5*cos(2.0*pi*iTime/10.0))*pi/6.0;
    float c = cos(mu);
    float s = sin(mu);
    mat2 M = mat2(c,s,-s,c);

    //if(c*(p.y*c-p.x*s)<0.0) return 0.0; // through line
    //if(abs(c*(p.y*c-p.x*s))>s) return 0.0; // envelope
    //if(abs((p.y*c/s))>n+1.0) return 0.0; // bound-y
    //if(abs(p.x)>n+1.0) return 0.0; // bound-x
    //if(box(p,vec2(n+1.0,(n+1.0)/c*s))>0.0) return 0.0; // bounding box

    p = p*M;

    p *= c;

    p.x -= clamp(round(p.x),-n,+n); // from iq's https://www.shadertoy.com/view/3syGzz limited rpetition SDF
    vec2 p0 = vec2(p.x-0.5,p.y);
    vec2 p1 = vec2(p.x+0.5,p.y);

    //return min(box(M*p0-vec2(0,s),0.5*vec2(c,0)),box(M*p1-vec2(0,s),0.5*vec2(c,0)))/c; // line variant
    return min(box(M*p0,0.5*vec2(c,s)),box(M*p1,0.5*vec2(c,s)))/c; // box variant
    return // triangular variant
        min
        (
            hseg(p),
            min
            (
                max(-p.y,box(M*p0,0.4*vec2(c,s))),
                max(-p.y,box(M*p1,0.4*vec2(c,s)))
            )
        )
        /
        c;
}

float map(vec2 p)
{
    return adaptive_stairs(p, 2.0);
}

void mainImage(out vec4 o, vec2 I)
{
    vec2 R = iResolution.xy;

    vec2 p = 2.5*(I+I-R)/R.y;
    vec2 m = 2.5*(2.0*iMouse.xy-R)/R.y;

    float d = map(p);

    // iq's coloring
    vec3 col = (d<0.0) ? vec3(0.6,0.8,1.0) : vec3(0.9,0.6,0.3);
    col *= 1.0 - exp(-9.0*abs(d));
    col *= 1.0 + 0.2*cos(128.0*abs(d));
    col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.015,abs(d)) );
    // iq's mouse distance visualization
    if( iMouse.z>0.001 )
    {
        d = map(m);
        col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
        col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

    o = vec4(col, 1.0);
}
