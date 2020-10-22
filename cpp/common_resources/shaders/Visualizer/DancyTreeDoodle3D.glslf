

//https://www.shadertoy.com/view/4lVyzh

#version 330 core


in       vec2      fUV;
uniform sampler2D  iChannel0; 
out vec4 gl_FragColor;

uniform float iTime;
//uniform float dt;
uniform vec2  iResolution;



#define pi 3.1415926
vec3 light;

float ln (vec3 p, vec3 a, vec3 b, float R) { 
    float r = dot(p-a,b-a)/dot(b-a,b-a);
    r = clamp(r,0.,1.);
    p.x+= 0.2*sqrt(R)*smoothstep(1.,0.,abs(r*2.-1.))*cos(pi*(2.*iTime));
    return length(p-a-(b-a)*r)-R*(1.5-0.4*r);
}

mat2 ro (float a) {
	float s = sin(a), c = cos(a);
    return mat2(c,-s,s,c);
}

float map (vec3 p) {
    float l = length(p-light)-1e-2;
    l = min(l,abs(p.y+0.4)-1e-2);
    l = min(l,abs(p.z-0.4)-1e-2);
    l = min(l,abs(p.x-0.7)-1e-2);
    p.y += 0.4;
    p.z += 0.1;
    p.zx *= ro(.5*iTime);
    vec2 rl = vec2(0.02,.25+ 0.01*sin(pi*4.*iTime));
    for (int i = 1; i < 11; i++) {
        
        l = min(l,ln(p,vec3(0),vec3(0,rl.y,0),rl.x));
    	p.y -= rl.y;
        p.xy *= ro(0.2*sin(3.1*iTime+float(i))+sin(0.222*iTime)*(-0.1*sin(0.4*pi*iTime)+sin(0.543*iTime)/max(float(i),2.)));
        p.x = abs(p.x);
        p.xy *= ro(0.6+0.4*sin(iTime)*sin(0.871*iTime)+0.05*float(i)*sin(2.*iTime));
        p.zx *= ro(0.5*pi+0.2*sin(0.5278*iTime)+0.8*float(i)*(sin(0.1*iTime)*(sin(0.1*pi*iTime)+sin(0.333*iTime)+0.2*sin(1.292*iTime))));
        
        rl *= (.7+0.015*float(i)*(sin(iTime)+0.1*sin(4.*pi*iTime)));
        
        l=min(l,length(p)-0.15*sqrt(rl.x));
    }
	return l;
}

vec3 march (vec3 p, vec3 d) {
    float o = 1e3;
    for (int i = 0; i < 24; i++) {
        float l = map(p);
    	p += l*d;
        if (l < 1e-3)break;
    }
    return p;
}

vec3 norm (vec3 p) { // iq
    vec2 e = vec2 (.001,0.);
    return normalize(vec3(
            map(p+e.xyy) - map(p-e.xyy),
            map(p+e.yxy) - map(p-e.yxy),
            map(p+e.yyx) - map(p-e.yyx)
        ));
}

//void mainImage( out vec4 C, in vec2 U )
void main( ){

    vec4 C;
    light = vec3(0.2*sin(iTime),0.5,-.5);
    //if (iMouse.z > 0.) light = vec3(vec2(-0.5,0.5)*0.+0.7*(iMouse.xy-0.5*R)/R.y,-.3);
    
    vec2 U = (fUV-0.5)*2.0;
    vec2 R = iResolution.xy;
    U.x*=(iResolution.x/iResolution.y);
    //U = (U-0.5*R)/R.y;

    vec3 p = vec3(0,0,-1);
    vec3 d = normalize(vec3(U,1));
    p =  march(p,d);
    vec3 n = norm(p);
	C = 0.6+0.4*sin(1.1*vec4(1,2,3,4)*dot(d,n));
    vec3 D = light-p;
    d = normalize(D);
    vec3 lp = march(p+d*1e-2,d);
    C *= 2.5*(dot(d,n))*(.3+0.7*length(lp-p)/length(light-p));
    C = atan(C)/pi*2.;

    gl_FragColor = C;
}


