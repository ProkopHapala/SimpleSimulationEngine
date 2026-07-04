
## Common

```glsl
float getOccupancy(vec2 uv) { return texture(iChannel0, uv).r; }

bool isIn(vec2 uv, float threshold) { return getOccupancy(uv) > threshold; }

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // This buffer averages out multiple signed distance fields
    // wherein the threshold that deliniates inside from outside
    // varies (from 0 to 1) between different signed distance
    // computations.
    // This essentially allows a signed distance field
    // where there is no binary distinction between inside and outside,
    // where "low occupancy" pixels have lower influence
    // rather than having no influence at all.

    // Compute the "noise" which is used to generate
    // and combine multiple samples.
    // The noise vector contains:
    //     x, y: offset to sampling position (antialiasing)
    //     z:    threshold deliniating inside from outside,
    //           as the input texture contains a range of values
    //           from 0 to 1 and must be converted to binary
    //           at a given threshold.
    vec3 noise;
    {
        noise = vec3(0);
        
        // Calculate oising over time.
        vec3 temporalNoise;
        {
            // Use the golden ratio as it should land
            // on all fractional values eventually.
            temporalNoise = vec3(iFrame, iFrame+1, iFrame+2);
            temporalNoise *= 1.618033;
            
            // We floor this one early to prevent
            // loss of precision when iFrame becomes large.
            temporalNoise -= floor(noise);
        }
        noise += temporalNoise;
    
        #ifdef SPATIAL_NOISE
        // Add noising over space.
        // (Currently disabled; messes up the
        // gradient and especially the curvature
        // of the resulting map.)
        vec3 spatialNoise;
        {
            // Noise is added to vary the threshold
            // per pixel to speed up apparent convergence,
            // but the converged result shouldn't change.
            
            vec2 noiseUV = fragCoord.xy / iChannelResolution[1].xy;
            spatialNoise = texture(iChannel1, noiseUV).r;
        }
        noise += spatialNoise;
        #endif
        
        // Wrap values around from 0 to 1.
        noise -= floor(noise);
        
        // Center the sampling position offset
        // to a range within -0.5 to +0.5.
        noise.xy -= 0.5f;
    }
    
    // Compute the signed distance to an edge
    // (limited to a certain search distance)
    // from -1 to 1.
    float distanceToEdge;
    {
        // Sample with "jitter" for accumulative antialiasing.
        vec2 samplingCenter = fragCoord + noise.xy;
        vec2 samplingCenterUV = samplingCenter / iChannelResolution[0].xy;

        const int iRange = 16;
        const float range = float(iRange);
        const float maxSqrDist = range*range;
        vec2 startPosition = samplingCenter;

        // We need to know whether our center pixel,
        // the one we're currently calculating the signed distance for,
        // is inside the volume or outside,
        // as this determines its sign:
        // Inside -> negative
        // Outside -> positive
        bool fragIsIn = isIn(samplingCenterUV, noise.z);
        
        float squaredDistanceToEdge = maxSqrDist;
        for(int dx=-iRange; dx <= iRange; dx++)
        {
            for(int dy=-iRange; dy <= iRange; dy++)
            {
                vec2 delta = vec2(dx, dy);
                vec2 scanPosition = startPosition + vec2(dx, dy);
                float scanDistanceSqr = dot(delta, delta);

                // Ideally we'd use a precomputed sampling pattern
                // that avoids testing the corners that are out-of-range.
                // Perhaps the compiler already unrolls the for loop 
                // and culls them?
                if(scanDistanceSqr >= maxSqrDist)
                    continue;
                
                // Already found one closer? Skip.
                if(scanDistanceSqr >= squaredDistanceToEdge)
                    continue;

                bool scanIsIn = isIn(scanPosition / iChannelResolution[0].xy, noise.z);
                
                // Sign change?
                if (scanIsIn != fragIsIn)
                {
                    // We found a boundary!
                    squaredDistanceToEdge = scanDistanceSqr;
                }
            }
        }

        distanceToEdge = sqrt(squaredDistanceToEdge);

        // The minimum distance is always 1,
        // but the boundary lies between the pixels halfway.
        // Correct the discontinuity this creates
        // in the distance field (a flat region at the sign change
        // boundary) by subtracting the pixel radius
        // from the distance.
        distanceToEdge -= 0.5;

        // Make the distance signed:
        // Inside -> negative
        // Outside -> positive
        distanceToEdge = fragIsIn ? -distanceToEdge : distanceToEdge;
        
        // Convert distance in pixels 
        // from -range to +range
        // to 0 to 1 
        distanceToEdge /= range * 2.;
        distanceToEdge = 0.5 - distanceToEdge;
    }
    
    // This step will prevent us from generating
    // a true signed distance field, but is useful
    // for terrain because it removes discontinuities
    // in the gradient.
    distanceToEdge = smoothstep(0., 1., distanceToEdge);
    
    // Optionally repeating this step will remove
    // curvature discontinuities, but
    // reduces the influence of distant pixels too much.
    //distanceToEdge = smoothstep(0., 1., distanceToEdge);
    
    // NOTE: you can apply other transformations to
    // "distanceToEdge" to change the character of the terrain.
    // For example, squaring the values makes peaks pointier
    // and cracks flatter.
    // distanceToEdge *= distanceToEdge;
    
    fragColor = vec4(distanceToEdge, distanceToEdge, distanceToEdge, 1.0);

    // Accumulate samples over time.
    {
        vec2 uv = (fragCoord.xy) / iChannelResolution[2].xy;
        vec4 oldColor = texture(iChannel2, uv);

        // This block of code was to make the terrain
        // "grow" instead of expose the poorly resolved
        // initial state with few samples.
        /*if(oldColor.a == 0.)
        {
            oldColor.rgb = vec3(25.);
            oldColor.a = 50.;
        }*/

        fragColor += oldColor;
    }
}

```

## Buffer A

```glsl
float perlin(vec2 uv)
{
    uv += vec2(8.2813, 1.42114);
    uv /= 4.0;
    vec2 occ = vec2(0);
    float a = 1.0;
    for(int i = 0; i < 7; i++)
    {
        occ += vec2(texture(iChannel0, uv).r, 1) * a;
        uv *= 0.5;
        a *= 2.0;
    }
    float v = occ.x / occ.y;
    
    // Increase contrast
    // (though this will flatten the tops and bottoms).
    v = v * 2.0 - 1.0;
    v = tanh(v * 2.0);
    v = v * 0.5 + 0.5;
    
    // Make peaks pointer and valleys flatter.
    v *= v;
    
    return v;
       
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = fragCoord.xy / iChannelResolution[0].xy;
    float height = perlin(uv);
    fragColor = vec4(height, height, height, 1);
}
```

## Image

```glsl
vec2 compress(vec2 vec)
{
    float mag = sqrt(dot(vec, vec));

    float newMag = mag;

    // Softly compress the range 
    // from 0 to +Inf
    // to 0 to +1 
    // instead of clipping values
    // when contrast is used for visualization.
    newMag = newMag / (0.5 + newMag);

    vec *= newMag / (mag + 0.00001);
    return vec;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uvC = fragCoord.xy / iResolution.xy;
	vec2 uvR = (fragCoord.xy + vec2(1,0)) / iResolution.xy;
	vec2 uvT = (fragCoord.xy + vec2(0,1)) / iResolution.xy;
	vec2 uvL = (fragCoord.xy - vec2(1,0)) / iResolution.xy;
	vec2 uvB = (fragCoord.xy - vec2(0,1)) / iResolution.xy;
    vec2 colC = texture(iChannel0, uvC).ra;
    vec2 colR = texture(iChannel0, uvR).ra;
    vec2 colT = texture(iChannel0, uvT).ra;
    vec2 colL = texture(iChannel0, uvL).ra;
    vec2 colB = texture(iChannel0, uvB).ra;
    
    colC.r /= colC.g;
    colR.r /= colR.g;
    colT.r /= colT.g;
    colL.r /= colR.g;
    colB.r /= colT.g;
    
    
    vec2 gradient = vec2(colL.r - colR.r, colB.r - colT.r);
    vec2 curvature = vec2
    (
        colL.r + colR.r - 2.0 * colC.r, 
        colB.r + colT.r - 2.0 * colC.r 
    ) * 0.5;
    
    
    // Enhance gradient+curvature contrast.
    vec2 both = compress(gradient * 50.0 - curvature * 1000.0);
    
    
    // Unmodified gradients:
    // vec2 gradient = vec2(colC.r - colR.r, colC.r - colT.r);
    
    fragColor = vec4(both * 0.5 + 0.5, colC.r, 1);
}

```
