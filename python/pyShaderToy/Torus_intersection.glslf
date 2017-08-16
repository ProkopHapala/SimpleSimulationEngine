// The MIT License
// Copyright © 2014 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



// Analytic intersection of a torus (degree 4 equation). Motivated by Antonalog's 
// shader (https://www.shadertoy.com/view/XdSGWy), and simplified the geometrically 
// impossible cases, and optimized coefficients. One can probably do better than 
// this though...
	

// f(x) = (|x|² + R² - r²)² - 4·R²·|xy|² = 0

float iTorus( in vec3 ro, in vec3 rd, in vec2 torus ){
	float Ra2 = torus.x*torus.x;
	float ra2 = torus.y*torus.y;
	float m = dot(ro,ro);
	float n = dot(ro,rd);
	float k  = (m - ra2 - Ra2)/2.0;
	float a  = n;
	float b  = n*n + Ra2*rd.z*rd.z + k;
	float c  = k*n + Ra2*ro.z*rd.z;
	float d  = k*k + Ra2*ro.z*ro.z - Ra2*ra2;
    //----------------------------------
    float a2 =  a*a;
	float p  =     -3.0*a2 + 2.0*b;
	float q  =  2.0*(a*( a2 - b) + c);
	float r  = a2*(-3.0*a2 + 4.0*b) + 4.0*(-2.0*a*c + d);
	p /= 3.0; r /= 3.0;
	float p2 = p*p;
	float Q  = p2 + r;
	float R  = (3.0*r - p2)*p - q*q;
	float h = R*R - Q*Q*Q;
	float z = 0.0;
	if( h < 0.0 ){
		float sQ = sqrt(Q);
		z = 2.0*sQ*cos( acos(R/(sQ*Q)) / 3.0 );
	}else{
		float sQ = pow( sqrt(h) + abs(R), 1.0/3.0 );
		z = sign(R)*abs( sQ + Q/sQ );
	}
	z = p - z;
    //----------------------------------
	float d1 = z   - 3.0*p;
	float d2 = z*z - 3.0*r;
	if( abs(d1)<1.0e-4 ){
		if( d2<0.0 ) return -1.0;
		d2 = sqrt(d2);
	}else{
		if( d1<0.0 ) return -1.0;
		d1 = sqrt( d1/2.0 );
		d2 = q/d1;
	}
    //----------------------------------
	float result = 1e20;
	h = d1*d1 - z + d2;
	if( h>0.0 ){
		h = sqrt(h);
		float t1 = -d1 - h - a;
		float t2 = -d1 + h - a;
		     if( t1>0.0 ) result=t1;
		else if( t2>0.0 ) result=t2;
	}
	h = d1*d1 - z - d2;
	if( h>0.0 ){
		h = sqrt(h);
		float t1 = d1 - h - a;
		float t2 = d1 + h - a;
		     if( t1>0.0 ) result=min(result,t1);
		else if( t2>0.0 ) result=min(result,t2);
	}
	return result;
}

// df(x)/dx
vec3 nTorus( in vec3 pos, vec2 tor )
{
	return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 p = (-iResolution.xy + 2.0*fragCoord.xy) / iResolution.y;

     // camera movement	
	float an = 0.5*iTime;
	vec3 ro = vec3( 2.5*cos(an), 1.0, 2.5*sin(an) );
    vec3 ta = vec3( 0.0, 0.0, 0.0 );
    // camera matrix
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
	// create view ray
	vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

    // raytrace
	
	// raytrace-plane
	vec2 torus = vec2(1.0,0.5);
	float t = iTorus( ro, rd, torus );

    // shading/lighting	
	vec3 col = vec3(0.0);
	if( t>0.0 && t<100.0 ){
	    vec3 pos = ro + t*rd;
		vec3 nor = nTorus( pos, torus );
		float dif = clamp( dot(nor,vec3(0.57703)), 0.0, 1.0 );
		float amb = clamp( 0.5 + 0.5*dot(nor,vec3(0.0,1.0,0.0)), 0.0, 1.0 );
		col = vec3(0.2,0.3,0.4)*amb + vec3(1.0,0.9,0.7)*dif;
		col *= 0.8;
		
		col *= (0.8 + 0.2*(sin(pos*100.0)*0.5 + 0.5));
	}
	
	col = sqrt( col );
	
	fragColor = vec4( col, 1.0 );
}
