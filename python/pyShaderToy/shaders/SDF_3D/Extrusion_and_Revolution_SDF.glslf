// The MIT License
// Copyright © 2018 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


// Two techniques used to generate 3D primitives from 2D primitives - extrusion
// and Revolution. Assuming the 2D shape defines exact distances, the resulting
// 3D shape is exact and way cheaper to evaluate than producing the same result
// with boolean operations.


// Related techniques:
//
// Elongation  : https://www.shadertoy.com/view/Ml3fWj
// Rounding    : https://www.shadertoy.com/view/Mt3BDj
// Onion       : https://www.shadertoy.com/view/MlcBDj
// Metric      : https://www.shadertoy.com/view/ltcfDj
// Combination : https://www.shadertoy.com/view/lt3BW2
// Repetition  : https://www.shadertoy.com/view/3syGzz
// Extrusion2D : https://www.shadertoy.com/view/4lyfzw
// Revolution2D: https://www.shadertoy.com/view/4lyfzw
//
// More information here: https://iquilezles.org/articles/distfunctions

//-------------------------------------------------

float opExtrusion( in vec3 p, in float sdf, in float h )
{
    vec2 w = vec2( sdf, abs(p.z) - h );
  	return min(max(w.x,w.y),0.0) + length(max(w,0.0));
}


vec2 opRevolution( in vec3 p, float w )
{
    return vec2( length(p.xz) - w, p.y );
}

//-------------------------------------------------



// https://iquilezles.org/articles/distfunctions2d
float sdCross( in vec2 p, in vec2 b, float r ) 
{
    p = abs(p); p = (p.y>p.x) ? p.yx : p.xy;
    
	vec2  q = p - b;
    float k = max(q.y,q.x);
    vec2  w = (k>0.0) ? q : vec2(b.y-p.x,-k);
    
    return sign(k)*length(max(w,0.0)) + r;
}

// https://iquilezles.org/articles/distfunctions2d
float sdVesica(vec2 p, float r, float d)
{
    p = abs(p);

    float b = sqrt(r*r-d*d); // can delay this sqrt
    return ((p.y-b)*d > p.x*b) 
            ? length(p-vec2(0.0,b))
            : length(p-vec2(-d,0.0))-r;
}

//---------------------------------

float map(in vec3 pos)
{
    float d = 1e10;
    

    // vesica 2D
    {
    vec3 q = pos - vec3(-3.0,0.0,-1.0);
    d = min(d,opExtrusion( q, sdVesica( q.xy, 0.7, 0.3 ), 0.005 ));
    }
    
    // extruded vesica
    {
    vec3 q = pos - vec3(-3.0,0.0,1.0);
    d = min(d,opExtrusion( q, sdVesica( q.xy, 0.7, 0.3 ), 0.4 ));
    }

    // cross 2D
    {
    vec3 q = pos - vec3(-1.0,0.0,-1.0);
    d = min(d,opExtrusion( q, sdCross( q.xy, vec2(0.8,0.35), 0.2 ), 0.005 ));
    }
    
    // extruded cross
    {
    vec3 q = pos - vec3(-1.0,0.0,1.0);
    d = min(d,opExtrusion( q, sdCross( q.xy, vec2(0.8,0.35), 0.2 ), 0.4 ));
    }

    // vesica 2D
    {
    vec3 q = pos - vec3(1.0,0.0,-1.0);
    d = min(d,opExtrusion( q, sdVesica( q.xy, 0.7, 0.3 ), 0.005 ));
    }
    
    // revolved vesica
    {
    vec3 q = pos - vec3(1.0,0.0,1.0);
     d = min( d, sdVesica(opRevolution(q,0.15-0.15*sin(iTime)), 0.7, 0.3 ) );
            
    }

    // cross 2D
    {
    vec3 q = pos - vec3(3.0,0.0,-1.0);
    d = min(d,opExtrusion( q, sdCross( q.xy, vec2(0.55,0.2), 0.1 ), 0.005 ));
    }
    
    // revolved cross
    {
    vec3 q = pos - vec3(3.0,0.0,1.0);
    d = min(d, sdCross( opRevolution(q,0.2+0.2*sin(iTime)), vec2(0.5,0.15), 0.1) );
    }

    return d;
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
    const float ep = 0.0001;
    vec2 e = vec2(1.0,-1.0)*0.5773;
    return normalize( e.xyy*map( pos + e.xyy*ep ) + 
					  e.yyx*map( pos + e.yyx*ep ) + 
					  e.yxy*map( pos + e.yxy*ep ) + 
					  e.xxx*map( pos + e.xxx*ep ) );
}

// https://iquilezles.org/articles/rmshadows
float calcSoftshadow( in vec3 ro, in vec3 rd, float tmin, float tmax, const float k )
{
	float res = 1.0;
    float t = tmin;
    for( int i=0; i<50; i++ )
    {
		float h = map( ro + rd*t );
        res = min( res, k*h/t );
        t += clamp( h, 0.02, 0.20 );
        if( res<0.005 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}


#define AA 2

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
   vec3 tot = vec3(0.0);
    
    #if AA>1
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        #else    
        vec2 p = (-iResolution.xy + 2.0*fragCoord)/iResolution.y;
        #endif
 
        vec3 ro = vec3(0.0,3.0,6.0);
        vec3 rd = normalize(vec3(p-vec2(0.0,1.0),-2.0));

        float t = 5.0;
        for( int i=0; i<64; i++ )
        {
            vec3 p = ro + t*rd;
            float h = map(p);
            if( abs(h)<0.001 || t>10.0 ) break;
            t += h;
        }

        vec3 col = vec3(0.0);

        if( t<10.0 )
        {
            vec3 pos = ro + t*rd;
            vec3 nor = calcNormal(pos);
            vec3  lig = normalize(vec3(1.0,0.8,0.2));
            float dif = clamp(dot(nor,lig),0.0,1.0);
            float sha = calcSoftshadow( pos, lig, 0.001, 1.0, 32.0 );
            col = vec3(0.025,0.05,0.08) + dif*sha*vec3(1.0,0.9,0.8);
        }

        col = sqrt( col );
	    tot += col;
    #if AA>1
    }
    tot /= float(AA*AA);
    #endif

	fragColor = vec4( tot, 1.0 );
}