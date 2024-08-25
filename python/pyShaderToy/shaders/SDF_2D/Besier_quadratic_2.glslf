// The MIT License
// Copyright © 2018 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Distance to a quadratic bezier segment, which can be solved analyically with a cubic.

// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and iquilezles.org/articles/distfunctions2d


// 0: exact, using a cubic colver
// 1: approximated
#define METHOD 0



float dot2( vec2 v ) { return dot(v,v); }
float cro( vec2 a, vec2 b ) { return a.x*b.y-a.y*b.x; }
float cos_acos_3( float x ) { x=sqrt(0.5+0.5*x); return x*(x*(x*(x*-0.008972+0.039071)-0.107074)+0.576975)+0.5; } // https://www.shadertoy.com/view/WltSD7

#if METHOD==0
// signed distance to a quadratic bezier
float sdBezier( in vec2 pos, in vec2 A, in vec2 B, in vec2 C, out vec2 outQ )
{    
    vec2 a = B - A;
    vec2 b = A - 2.0*B + C;
    vec2 c = a * 2.0;
    vec2 d = A - pos;

    // cubic to be solved (kx*=3 and ky*=3)
    float kk = 1.0/dot(b,b);
    float kx = kk * dot(a,b);
    float ky = kk * (2.0*dot(a,a)+dot(d,b))/3.0;
    float kz = kk * dot(d,a);      

    float res = 0.0;
    float sgn = 0.0;

    float p  = ky - kx*kx;
    float q  = kx*(2.0*kx*kx - 3.0*ky) + kz;
    float p3 = p*p*p;
    float q2 = q*q;
    float h  = q2 + 4.0*p3;


    if( h>=0.0 ) 
    {   // 1 root
        h = sqrt(h);
        vec2 x = (vec2(h,-h)-q)/2.0;

        #if 0
        // When p≈0 and p<0, h-q has catastrophic cancelation. So, we do
        // h=√(q²+4p³)=q·√(1+4p³/q²)=q·√(1+w) instead. Now we approximate
        // √ by a linear Taylor expansion into h≈q(1+½w) so that the q's
        // cancel each other in h-q. Expanding and simplifying further we
        // get x=vec2(p³/q,-p³/q-q). And using a second degree Taylor
        // expansion instead: x=vec2(k,-k-q) with k=(1-p³/q²)·p³/q
        if( abs(p)<0.001 )
        {
            float k = p3/q;              // linear approx
          //float k = (1.0-p3/q2)*p3/q;  // quadratic approx 
            x = vec2(k,-k-q);  
        }
        #endif

        vec2 uv = sign(x)*pow(abs(x), vec2(1.0/3.0));
        float t = uv.x + uv.y;

		// from NinjaKoala - single newton iteration to account for cancellation
        t -= (t*(t*t+3.0*p)+q)/(3.0*t*t+3.0*p);
        
        t = clamp( t-kx, 0.0, 1.0 );
        vec2  w = d+(c+b*t)*t;
        outQ = w + pos;
        res = dot2(w);
    	sgn = cro(c+2.0*b*t,w);
    }
    else 
    {   // 3 roots
        float z = sqrt(-p);
        #if 0
        float v = acos(q/(p*z*2.0))/3.0;
        float m = cos(v);
        float n = sin(v);
        #else
        float m = cos_acos_3( q/(p*z*2.0) );
        float n = sqrt(1.0-m*m);
        #endif
        n *= sqrt(3.0);
        vec3  t = clamp( vec3(m+m,-n-m,n-m)*z-kx, 0.0, 1.0 );
        vec2  qx=d+(c+b*t.x)*t.x; float dx=dot2(qx), sx=cro(a+b*t.x,qx);
        vec2  qy=d+(c+b*t.y)*t.y; float dy=dot2(qy), sy=cro(a+b*t.y,qy);
        if( dx<dy ) {res=dx;sgn=sx;outQ=qx+pos;} else {res=dy;sgn=sy;outQ=qy+pos;}
    }
    
    return sqrt( res )*sign(sgn);
}
#else

// This method provides just an approximation, and is only usable in
// the very close neighborhood of the curve. Taken and adapted from
// http://research.microsoft.com/en-us/um/people/hoppe/ravg.pdf
float sdBezier( vec2 p, vec2 v0, vec2 v1, vec2 v2, out vec2 outQ )
{
	vec2 i = v0 - v2;
    vec2 j = v2 - v1;
    vec2 k = v1 - v0;
    vec2 w = j-k;

	v0-= p; v1-= p; v2-= p;
    
	float x = cro(v0, v2);
    float y = cro(v1, v0);
    float z = cro(v2, v1);

	vec2 s = 2.0*(y*j+z*k)-x*i;

    float r =  (y*z-x*x*0.25)/dot2(s);
    float t = clamp( (0.5*x+y+r*dot(s,w))/(x+y+z),0.0,1.0);
    
    vec2 d = v0+t*(k+k+t*w);
    outQ = d + p;
	return length(d);
}
#endif

float udSegment( in vec2 p, in vec2 a, in vec2 b )
{
	vec2 pa = p - a;
	vec2 ba = b - a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 m = (2.0*iMouse.xy-iResolution.xy)/iResolution.y;
    
	vec2 v0 = vec2(1.3,0.9)*cos(iTime*0.5 + vec2(0.0,5.0) );
    vec2 v1 = vec2(1.3,0.9)*cos(iTime*0.6 + vec2(3.0,4.0) );
    vec2 v2 = vec2(1.3,0.9)*cos(iTime*0.7 + vec2(2.0,0.0) );
    
    vec2 kk;
    
    float d = sdBezier( p, v0,v1,v2, kk ); 
    
    float f = smoothstep(-0.2,0.2,cos(2.0*iTime));
    vec3 col = (d>0.0) ? vec3(0.9,0.6,0.3) : vec3(0.65,0.85,1.0);
	col *= 1.0 - exp(-4.0*abs(d));
	col *= 0.8 + 0.2*cos(110.0*d);
	col = mix( col, vec3(1.0), 1.0-smoothstep(0.0,0.015,abs(d)) );
    
    if( iMouse.z>0.001 )
    {
    vec2 q;
    d = sdBezier(m, v0,v1,v2, q ); 
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, abs(length(p-m)-abs(d))-0.0025));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-m)-0.015));
    col = mix(col, vec3(1.0,1.0,0.0), 1.0-smoothstep(0.0, 0.005, length(p-q)-0.015));
    }
 
    if( cos(0.5*iTime)<-0.5 )
    {
        d = min( udSegment(p,v0,v1),
                 udSegment(p,v1,v2) );
        d = min( d, length(p-v0)-0.02 );
        d = min( d, length(p-v1)-0.02 );
        d = min( d, length(p-v2)-0.02 );
        col = mix( col, vec3(1,0,0), 1.0-smoothstep(0.0,0.007,d) );
    }

	fragColor = vec4(col,1.0);
}