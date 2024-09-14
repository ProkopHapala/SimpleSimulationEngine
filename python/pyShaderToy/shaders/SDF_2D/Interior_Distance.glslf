// The MIT License
// Copyright © 2020 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// This shader shows how the usual union SDF operator implemented
// with min() produces the incorrect SDF in the interior of the shapes.
// If the modeling is reversed and the exterior of the negative
// space of the shape is modeled, then the exterior distance is wrong.
//
// While this is not an issue in most cases in practice, it can be
// a problem in shaders that need to raymarch or do some other
// distance based volumetric operations in the interior of the
// shapes.
//
// This shader implements both correct interior and exterior distances
// by modeling the boundary of the shape instead. Alternativelly, the
// appropriate interior or exterior correct SDF could be chosen as
// needed, but then there is double modeling work.
//
// More info here:
// https://iquilezles.org/articles/interiordistance



float dot2( in vec2 v ) { return dot(v,v); }
float msign( in float x ) { return (x>0.0)?1.0:-1.0; }

// https://iquilezles.org/articles/distfunctions2d
float sdCircle( in vec2 p, in vec2 c, in float r )
{
    return length(p-c) - r;
}

// https://iquilezles.org/articles/distfunctions2d
float sdBox( in vec2 p, in vec2 c, in vec2 b ) 
{
    vec2 q = abs(p-c) - b;
    return min(max(q.x,q.y),0.0) + length(max(q,0.0));
}

// https://iquilezles.org/articles/distfunctions2d
vec2 sdSqLine( in vec2 p, in vec2 a, in vec2 b )
{
	vec2 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return vec2( dot2(pa-ba*h), ba.x*pa.y-ba.y*pa.x );
}

float sdCrescent(vec2 p, float r0, float r1, float d, float sign0, float sign1)
{
    float a = (r0*r0 - r1*r1 + d*d) / (2.0 * d);
    
    if( a < r0)
    {
        p.y = abs(p.y);
        float b = sqrt(r0*r0-a*a);
        float k = p.y*a - p.x*b;
        float h = min(d*sign0*(d*(p.y-b)-k ),
                      d*sign1*k);
        if (h>0.0)
        {
            return length(p-vec2(a,b));
        }
    }
    
    return max(sign0*(length(p          )-r0),
               sign1*(length(p-vec2(d,0))-r1));
}

// https://iquilezles.org/articles/distfunctions2d
vec2 sdSqArc( in vec2 p, in vec2 a, in vec2 b, in float h, float d2min )
{
    vec2  ba  = b-a;
    float l   = length(ba);
    float ra2 = h*h + l*l*0.25;

    // recenter
    p -= (a+b)/2.0 + vec2(-ba.y,ba.x)*h/l;
    
    float m = ba.y*p.x-ba.x*p.y;
    float n = dot(p,p);
    
    if( abs(h)*abs(ba.x*p.x+ba.y*p.y) < msign(h)*l*0.5*m )
    {
        d2min = min( d2min, n + ra2 - 2.0*sqrt(n*ra2) );
    }

    return vec2(d2min, -max(m,ra2-n) );
}


//------------------------------------------------------------



// SDF of a shape made of a set line and arc segments
float sdShape( in vec2 p, int kType[7], float kPath[17] )
{
    vec2 vb = vec2(kPath[0],kPath[1]);
    
    float d = dot2(p-vb);
    int off = 0;
    float s = 1.0;
    for( int i=0; i<kType.length(); i++ )
    {
        vec2 va = vb;
        vec2 ds;
        
        if( kType[i]==0) // line (x,y)
        {
            vb = vec2(kPath[off+2],kPath[off+3]);
            ds = sdSqLine( p, va, vb );
            off += 2;
        }
        else if( kType[i]==1) // arc (x,y,r)
        {
            vb = vec2(kPath[off+3],kPath[off+4]);
            ds = sdSqArc(p, va, vb, kPath[off+2], d );
        	off += 3;

        }
        
        // in/out test
        bvec3 cond = bvec3( p.y>=va.y, p.y<vb.y, ds.y>0.0 );
        if( all(cond) || all(not(cond)) ) s*=-1.0;  

        d = min( d, ds.x );
    }
	return s*sqrt(d);
}
              
// correct outside, incorrect inside
float sdA( in vec2 p )
{
    float d = sdCircle( p, vec2(-0.4, 0.3), 0.5);
    d = min(d,sdBox( p, vec2( 0.4,-0.4), vec2(0.4,0.4) ));
    d = min(d,sdBox( p, vec2( 0.0, 0.0), vec2(0.4,0.8) ));
    return d;
}

// correct inside, incorrect outside
float sdB( in vec2 p )
{
   float d =     sdBox( p, vec2( 0.0, 1.0), vec2(2.0,0.2) );
       d = min(d,sdBox( p, vec2( 1.2, 1.0), vec2(0.8,1.0) ));
       d = min(d,sdBox( p, vec2( 1.4,-0.3), vec2(0.6,0.9) ));
       d = min(d,sdBox( p, vec2( 0.0,-1.0), vec2(1.0,0.2) ));
       d = min(d,sdBox( p, vec2(-1.2,-0.8), vec2(0.8,0.6) ));
       d = min(d,sdBox( p, vec2(-1.5, 0.3), vec2(0.6,0.7) ));
       d = min(d,sdCrescent( p-vec2(-0.4-1.0, 0.3), 1.1, 0.5, 1.0, 1.0, -1.0 ));
    return -d;
}

// correct both in side and outside
float sdC( in vec2 p )
{
	int kType[] = int[](0,0,0,0,0,0,1);
	float kPath[] = float[](-0.4, 0.8,
                             0.4, 0.8,
                             0.4,-0.0,
                             0.8,-0.0,
                             0.8,-0.8,
                            -0.4,-0.8,
                            -0.4,-0.2, 0.0,
                            -0.4, 0.8 );
    return sdShape(p,kType,kPath );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // normalized pixel coordinates
    vec2 p = (fragCoord*2.0-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
    
    // distance computations
    float dWrongInterior = sdA(p); // interior modeling
    float dWrongExterior = sdB(p); // exterior modeling
    float dCorrectBoth   = sdC(p); // boundary modeling

    // animation
    float f = fract(iTime/8.0);
    float g = fract(iTime/2.0);
    float d = (f<0.5) ? ((g<0.5)?dWrongInterior:dCorrectBoth) 
                      : ((g<0.5)?dWrongExterior:dCorrectBoth);
    
    // coloring
    vec3 col = (d<0.0) ? vec3(0.6,0.8,1.0) : vec3(0.9,0.6,0.3);
    col *= 1.0 - exp(-9.0*abs(d));
	col *= 1.0 + 0.2*cos(128.0*abs(d));
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.015,abs(d)) );

     // interactivity
    if( iMouse.z>0.001 )
    {
        float dWrongInterior = sdA(m); // interior modeling
        float dWrongExterior = sdB(m); // exterior modeling
        float dCorrectBoth   = sdC(m); // boundary modeling
        float d = (f<0.5) ? ((g<0.5)?dWrongInterior:dCorrectBoth) 
                          : ((g<0.5)?dWrongExterior:dCorrectBoth);

        col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
        col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    }

    // output
	fragColor = vec4(col, 1.0);
}