#version 330 core
// Generate a signed distance field (SDF) into a small RGBA32F texture
// Uses fullscreen quad VS (fs_quad.glslv). Outputs SDF in R channel (in texel units)
// and duplicates into other channels for easy inspection.

out vec4 color;
in vec2 uv;               // from fs_quad.glslv
uniform vec3 iResolution; // (W,H,1) of the target FBO/texture
uniform int  iFrame;      // not used here but provided for consistency
// Optional driver: (cx, cy, z, radius) in normalized [0,1] coords
uniform vec4 driver;

// Signed distance to a circle centered at c with radius r (in normalized [0,1] space)
float sd_circle(vec2 p, vec2 c, float r){ return length(p - c) - r; }

void main(){
    vec2 p = uv;                                      // [0,1]
    vec2 c = clamp(driver.xy, vec2(0.0), vec2(1.0));  // center in [0,1]
    float r = driver.w > 0.0 ? driver.w : 0.35;       // radius in [0,1]

    // --- GEN STAGE A: write UV to check FBO path ---------------------------
    // color = vec4(p, 0.0, 1.0); return;  // enable to verify FBO render works

    // --- GEN STAGE B: binary circle mask for visibility -------------------
    // float inside = length(p - c) < r ? 1.0 : 0.0; color = vec4(inside,inside,inside,1.0); return;

    // --- GEN STAGE C: normalized SDF preview (0..1) -----------------------
    // float sd_norm_dbg = sd_circle(p, c, r);
    // float d01 = clamp(sd_norm_dbg*8.0 + 0.5, 0.0, 1.0); // visualize ~±0.125 band
    // color = vec4(d01, d01, d01, 1.0); return;

    // Default: True SDF in texel units; robust defaults ensure visibility
    float sd_norm = sd_circle(p, c, r);               // normalized units
    float px = 1.0 / max(1.0, min(iResolution.x, iResolution.y));
    float sd_tex = sd_norm / px;                      // signed distance [texels]
    color = vec4(sd_tex, sd_tex, sd_tex, 1.0);        // store SDF in all channels
}
