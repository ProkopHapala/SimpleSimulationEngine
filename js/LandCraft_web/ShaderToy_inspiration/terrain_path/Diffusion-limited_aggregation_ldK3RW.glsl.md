
## Buffer A

```glsl
// 2016 David A Roberts <https://davidar.io>

// 1 out, 3 in... <https://www.shadertoy.com/view/4djSRW>
#define MOD3 vec3(.1031,.11369,.13787)
float hash13(vec3 p3) {
	p3 = fract(p3 * MOD3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.x + p3.y)*p3.z);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    fragColor = vec4(0,0,0,1);
    if(iFrame < 10 || iMouse.z > 0.) {
        if(length(fragCoord-10.) < 2.)
            fragColor.y = 1.;
        //else if(hash13(vec3(fragCoord,iFrame)) < 0.05)
        //    fragColor.x = 1.;
        return;
    }
    
    if(fragCoord.x < 1. || fragCoord.x > iResolution.x-1. ||
       fragCoord.y < 1. || fragCoord.y > iResolution.y-1.) {
        fragColor.x = 1.;
        return;
    }
    
    vec4 c  = texture(iChannel0, (fragCoord + vec2( 0, 0)) / iResolution.xy);
    vec4 n  = texture(iChannel0, (fragCoord + vec2( 0, 1)) / iResolution.xy);
    vec4 ne = texture(iChannel0, (fragCoord + vec2( 1, 1)) / iResolution.xy);
    vec4 e  = texture(iChannel0, (fragCoord + vec2( 1, 0)) / iResolution.xy);
    vec4 se = texture(iChannel0, (fragCoord + vec2( 1,-1)) / iResolution.xy);
    vec4 s  = texture(iChannel0, (fragCoord + vec2( 0,-1)) / iResolution.xy);
    vec4 sw = texture(iChannel0, (fragCoord + vec2(-1,-1)) / iResolution.xy);
    vec4 w  = texture(iChannel0, (fragCoord + vec2(-1, 0)) / iResolution.xy);
    vec4 nw = texture(iChannel0, (fragCoord + vec2(-1, 1)) / iResolution.xy);
    
    // aggregation
    fragColor.y = clamp(
        c.y + c.x * (n.y + ne.y + e.y + se.y + s.y + sw.y + w.y + nw.y), 0., 1.);
    
    bool nc = int(4.*hash13(vec3(fragCoord + vec2( 0, 1), iFrame))) == 0;
    bool ec = int(4.*hash13(vec3(fragCoord + vec2( 1, 0), iFrame))) == 1;
    bool sc = int(4.*hash13(vec3(fragCoord + vec2( 0,-1), iFrame))) == 2;
    bool wc = int(4.*hash13(vec3(fragCoord + vec2(-1, 0), iFrame))) == 3;
    
    // diffusion
    fragColor.x = clamp(
    	n.x * float(nc) + e.x * float(ec) + s.x * float(sc) + w.x * float(wc) +
    	floor((1. + (n.x + e.x + s.x + w.x)/100.) * hash13(vec3(fragCoord,iFrame))), 0., 1.);
    
    fragColor.x *= 1. - fragColor.y;
}
```


## Image

```glsl
// 2016 David A Roberts <https://davidar.io>

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = fragCoord.xy / iResolution.xy;
	fragColor = (texture(iChannel0,uv) + texture(iChannel1,uv) + texture(iChannel2,uv) + texture(iChannel3,uv)) / 5.;
}
```
