#define exact

//circle centers can be anywhere in this
float crescent(vec2 p, vec2 p0, vec2 p1, float r0, float r1, float sign0, float sign1) {
    
    p -= p0;
    p1 -= p0;
#ifdef exact
    float d = length(p1);
    vec2 n = normalize(p1);
    vec2 nt = vec2(n.y,-n.x);
    p += -nt*max(dot(p,nt),0.0)*2.0;
    
    float a = (r0*r0 - r1*r1 + d*d) / (2.0 * d);
    
    if (a < r0) {
        float b = sqrt(r0*r0-a*a);
        vec2 corner = a*n-b*nt;
        vec2 q = p-corner;

        vec2 n0 = d*(corner-p1);
        n0 = vec2(n0.y,-n0.x);

        vec2 n1 = d*corner;
        n1 = vec2(-n1.y,n1.x);

        if (min(sign0*dot(q,n0),sign1*dot(q,n1)) > 0.0) {
            return length(q);
        }
	}
#endif
    float len = sign0*(length(p)-r0);
    len = max(len,sign1*(length(p-p1)-r1));
    return len;
}

//simplified version by iq
float crescent2(vec2 p, float r0, float r1, float d, float sign0, float sign1)
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

float crescent(vec2 p, float r0, float r1, float d, float sign0, float sign1) {
    
    //taken from:
    //https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    float a = (r0*r0 - r1*r1 + d*d) / (2.0 * d);
    
#ifdef exact
    if (a < r0) {
    	p.y = abs(p.y);
        
        float b = sqrt(r0*r0-a*a);
        vec2 corner = vec2(a,b);
        vec2 q = p-corner;
        
        vec2 n0 = d*(corner-vec2(d,0));
        n0 = vec2(n0.y,-n0.x);
        
        vec2 n1 = d*corner;
        n1 = vec2(-n1.y,n1.x);
		
        //uses dot product to determine which side of the lines the pixel is in
        if (min(sign0*dot(q,n0),sign1*dot(q,n1)) > 0.0) {
            return length(q);
        }
        
    }
#endif
    
    return max(sign0*(length(p)-r0),sign1*(length(p-vec2(d,0))-r1));
}

float line( vec2 pa, vec2 ba) 
{
    float h = max( dot(pa,ba)/dot(ba,ba), 0.0 );
    return length( pa - ba*h );
}

float visualize(vec2 p, float r0, float r1, float d, float sign0, float sign1) {
    
    //taken from:
    //https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    float a = (r0*r0 - r1*r1 + d*d) / (2.0 * d);
    
    float len = abs(length(p)-r0);
    len = min(len,abs(length(p-vec2(d,0))-r1));
    
    if (a < r0) {
    	p.y = abs(p.y);
        float b = sqrt(r0*r0-a*a);
        vec2 corner = vec2(a,b);
        vec2 q = p-corner;

        vec2 n0 = normalize((corner-vec2(d,0)));
        //n0 = vec2(n0.y,-n0.x);

        vec2 n1 = normalize(corner);
        //n1 = vec2(-n1.y,n1.x);

        len = min(len,abs(line(q,sign1*n0*10.0)));
        len = min(len,abs(line(q,sign0*n1*10.0)));
    }
    
    return max(0.5-len*iResolution.y/16.0,0.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (fragCoord*2.0-iResolution.xy)/iResolution.y;
    
    float t = iTime;
    float r0 = 0.8;
    float r1 = 0.6+sin(t*0.8)*0.3;
    float d = sin(t)*1.2;
    
    float sign0 = sign(fract(iTime*0.1)-0.5);
    float sign1 = sign(fract(iTime*0.05)-0.5);
    
    //finding crescent shape distance
    float len = crescent(uv, r0, r1, d, sign0, sign1);
    
    //distance field coloring by iq https://www.shadertoy.com/view/4lcBWn
    //this is a bit different, i wanted the border to be resolution independent
    vec3 col = vec3(1.0) - sign(len)*vec3(0.1,0.4,0.7);
	col *= 1.0 - exp(-iResolution.y*0.008*abs(len));
	col *= 0.8 + 0.2*cos(iResolution.y*0.3*abs(len));
	col = mix( col, vec3(1.0), max(1.0-abs(len)*iResolution.y*0.2,0.0));
	
	fragColor = vec4(col*col, 1.0);
    
    fragColor += visualize(uv, r0, r1, d, sign0, sign1);
    
    //quadtree part
    /**
    vec2 p = uv;
    vec2 fp = floor(p);
    vec2 lp = p-fp;
    float size = 1.0;
    
    for (int i = 0; i < 5; i++) {
    	float len = crescent(fp+size*0.5, r0, r1, d, sign0, sign1);
        if (abs(len) > size*0.5*sqrt(2.0)) break;
        
        vec2 q = step(0.5,lp);
        lp = lp*2.0-q;
        size *= 0.5;
        fp += q*size;
        
    }
    lp = abs(lp-0.5);
    float a = max(1.0-(1.0-max(lp.x,lp.y)*2.0)*size*iResolution.y*0.125,0.0)*0.1;
    fragColor += vec4(a);
	/**/
    
    // squareroot for 2.0 gamma
    fragColor = sqrt(fragColor);
}