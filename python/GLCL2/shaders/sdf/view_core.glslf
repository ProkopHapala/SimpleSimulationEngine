#version 330 core
// View a signed distance field (SDF) stored in RGBA32F texture
// The SDF is in texel units in the R channel (others duplicated)

in vec2 uv;
out vec4 color;
uniform sampler2D iChannel0;  // bound to unit 0 by the framework
uniform vec3 iResolution;      // target framebuffer resolution (widget size)
uniform int  iFrame;
uniform vec4 driver;           // (cx, cy, zoom, unused)

void main(){
    // --- STAGE 1: UV sanity check -----------------------------------------
    // Uncomment to verify FS quad + UVs are correct (should show gradient)
    // color = vec4(uv, 0.0, 1.0);
    // return;

    // --- STAGE 2: Raw texture channel -------------------------------------
    // Uncomment to see raw red channel of sdfTex
    // float v = texture(iChannel0, uv).r;
    // color = vec4(v, v, v, 1.0);
    // return;

    // --- STAGE 3 (DEFAULT): SDF, no zoom, fixed AA ------------------------
    // This should clearly show a crisp circle regardless of driver.z
    //float d = texture(iChannel0, uv).r;           // signed distance [texels]
    //float w = 0.75;                                // fixed AA width [texels]
    //float a = smoothstep(-w, +w, -d);              // 1 inside (d<0), 0 outside
    //color = vec4(mix(vec3(0.0), vec3(1.0), a), 1.0);
    //color = vec4(mix(vec3(0.0), vec3(1.0), 1.-d), 1.0);


    //color = vec4( texture(iChannel0, uv).rgb, 1.0);
    vec3 rgb = texture(iChannel0, uv).rgb;
    color = vec4( rgb*rgb, 1.0);
    //color = vec4( sqrt(rgb), 1.0);
    //color = vec4( sqrt(sqrt(rgb)), 1.0);
    //color = vec4( floor((rgb*rgb*rgb)+0.5), 1.0);



    
    
    
    return;

    // --- STAGE 4: Zoom/pan SDF with robust AA -----------------------------
    // Enable by commenting out the 'return' above.
    // vec2 c = clamp(driver.xy, vec2(0.0), vec2(1.0));
    // float z = driver.z > 0.0 ? driver.z : 1.0;
    // vec2 suv = clamp((uv - c)/z + c, vec2(0.0), vec2(1.0));
    // float d2 = texture(iChannel0, suv).r;
    // float w2 = max(fwidth(d2), 0.75);
    // float a2 = smoothstep(-w2, +w2, -d2);
    // color = vec4(mix(vec3(0.0), vec3(1.0), a2), 1.0);
}
