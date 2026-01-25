
## Common

```glsl
struct quaternion
{
    vec4 value;
};

struct bounds
{
    vec3 mini, maxi;
};

bool rayIntersectBounds
(
    bounds b, 
    vec3 rayOrig, 
    vec3 rayDir,
    out vec3 nearHit,
    out vec3 farHit
)
{
    b.mini -= rayOrig;
    b.maxi -= rayOrig;
    
    vec3 signs = sign(rayDir);
    
   
    b.mini *= signs;
    b.maxi *= signs;
    rayDir *= signs;
    
    vec3 maxBounds = max((b.mini), (b.maxi));
    
    rayDir *= max(max(maxBounds.x, maxBounds.y),maxBounds.z) * 2.0; 
    
    
    
    {
        farHit = rayDir;
        vec3 clamped = min(farHit, maxBounds);
    
        vec3 scale = max(vec3(0.0), clamped / farHit);
        float minScale = min(min(scale.x, scale.y),scale.z); 
        farHit = (farHit * minScale) * signs + rayOrig;
    }
    
    /*{
        nearHit = -rayDir;
        vec3 clamped = clamp(nearHit, b.mini, b.maxi);
    
        vec3 scale = abs(clamped / nearHit);
        float minScale = max(max(scale.x, scale.y),scale.z); 
        nearHit = nearHit * minScale + rayOrig;
    }*/
    
    nearHit = rayOrig;

    return true;
}

quaternion FromAngleAxis(vec3 angleAxis)
{
    float mag = length(angleAxis);
    float halfAngle = mag * 0.5;
    float scalar = sin(halfAngle) / max(mag, 0.00001);
        
    quaternion q;
    q.value = vec4(angleAxis * scalar, cos(halfAngle));
    return q;
}

quaternion FromToRotation(vec3 from, vec3 to)
{
    vec3 xyz = cross(from, to);
    float w = sqrt(dot(from, from) * dot(to, to)) + dot(from, to);
    vec4 value = vec4(xyz, w);
    quaternion q;
    q.value = normalize(value);
    return q;
}

vec3 mul(quaternion q, vec3 v)
{
    vec3 t = 2.0 * cross(q.value.xyz, v);
    return v + vec3(q.value.w * t) + cross(q.value.xyz, t);
}

quaternion mul(quaternion a, quaternion b)
{
    quaternion q;
    q.value = vec4
    (
        a.value.wwww * b.value 
      + (a.value.xyzx * b.value.wwwx + a.value.yzxy * b.value.zxyy) * vec4(1.0, 1.0, 1.0, -1.0) 
      - a.value.zxyz * b.value.yzxz
    );
    return q;
}


//	Simplex 3D Noise 
//	by Ian McEwan, Ashima Arts
//
vec4 permute(vec4 x){return mod(((x*34.0)+1.0)*x, 289.0);}
vec4 taylorInvSqrt(vec4 r){return 1.79284291400159 - 0.85373472095314 * r;}

// Simplex 2D noise
//
vec3 permute(vec3 x) { return mod(((x*34.0)+1.0)*x, 289.0); }

float snoise(vec2 v){
  const vec4 C = vec4(0.211324865405187, 0.366025403784439,
           -0.577350269189626, 0.024390243902439);
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);
  vec2 i1;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;
  i = mod(i, 289.0);
  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
  + i.x + vec3(0.0, i1.x, 1.0 ));
  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy),
    dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;
  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;
  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 0.5+0.5*(130.0 * dot(m, g));
}

float snoise(float p){
	return snoise(vec2(p,0.));
}

float snoise(vec3 v){ 
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //  x0 = x0 - 0. + 0.0 * C 
  vec3 x1 = x0 - i1 + 1.0 * C.xxx;
  vec3 x2 = x0 - i2 + 2.0 * C.xxx;
  vec3 x3 = x0 - 1. + 3.0 * C.xxx;

// Permutations
  i = mod(i, 289.0 ); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients
// ( N*N points uniformly over a square, mapped onto an octahedron.)
  float n_ = 1.0/7.0; // N=7
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z *ns.z);  //  mod(p,N*N)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1), 
                                dot(p2,x2), dot(p3,x3) ) );
}

#define NUM_OCTAVES 16

float fbm(float x) {
	float v = 0.0;
	float a = 0.5;
	float shift = float(100);
	for (int i = 0; i < NUM_OCTAVES; ++i) {
		v += a * snoise(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}


float fbm(vec2 x) {
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(100);
	// Rotate to reduce axial bias
    mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	for (int i = 0; i < NUM_OCTAVES; ++i) {
		v += a * snoise(x);
		x = rot * x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}

float turbulent(vec2 x) {
	float v = fbm(x);
    return 1.-abs((v-0.5)*2.);
}



float fbm(vec3 x) {
	float v = 0.0;
	float a = 0.5;
	vec3 shift = vec3(100);
	for (int i = 0; i < NUM_OCTAVES; ++i) {
		v += a * snoise(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return (v+1.)/2.;
}

float turbulent(vec3 x) {
	float v = fbm(x);
    return 1.-abs((v-0.5)*2.);
}

float pattern( in vec2 p )
  {
    vec2 q = vec2( fbm( p + vec2(0.0,0.0) ),
                   fbm( p + vec2(5.2,1.3) ) );

    return fbm( p + 4.0*q );
  }

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

int radius = 1;
```

## Buffer A

```glsl
float getOccupancy(vec2 uv)
{
    return texture(iChannel0, uv).r;
}

bool isIn(vec2 uv, float threshold)
{
	return getOccupancy(uv) > threshold;
}

float squaredDistanceBetween(vec2 uv1, vec2 uv2)
{
    vec2 difference = uv1 - uv2;
    return dot(difference, difference);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = (fragCoord.xy) / iResolution.xy;
    
    //fragColor = vec4(vec3(getOccupancy(uv)), 1);
    //return;
    
    // Compute the noise.
    vec3 noise;
    {
        
        // Cycle the noise.
        // We do the cycling FIRST to prevent
        // loss of precision that might
        // round away the noise values.
        //   - Use the golden ratio as it should land
        //     on all fractional values eventually.
        noise = vec3(iFrame, iFrame+1, iFrame+2);
        noise *= 1.618033;
        noise -= floor(noise);
    
        // Noise is added to vary the threshold
        // per pixel to speed up apparent convergence,
        // but the converged result shouldn't change.
        // (Disable to prevent gradient noise).
        //vec2 noiseUV = fragCoord.xy / iChannelResolution[1].xy;
        //noise += texture(iChannel1, noiseUV).r;
        //noise -= floor(noise);
        
        noise.xy -= 0.5f;
    }
    
    // Compare with and without fuzziness.
    //if(uv.x > 0.5)
    //    noise.z = 0.5;
       
    vec2 samplingCenter = fragCoord + noise.xy;
    vec2 samplingCenterUV = samplingCenter / iResolution.xy;
    
    const float halfRange = 16.0;
    const int iRange = int(halfRange);
    const float maxSqrDist = halfRange*halfRange;
    vec2 startPosition = samplingCenter;
    
    bool fragIsIn = isIn(samplingCenterUV, noise.z);
    float squaredDistanceToEdge = maxSqrDist*2.;
    
    for(int dx=-iRange; dx <= iRange; dx++)
    {
        for(int dy=-iRange; dy <= iRange; dy++)
        {
            vec2 delta = vec2(dx, dy);
            vec2 circlizer = abs(delta);
            circlizer /= max(circlizer.x, circlizer.y);
            delta /= length(circlizer);
            vec2 scanPosition = startPosition + vec2(dx, dy);
            float scanDistance = squaredDistanceBetween(samplingCenter, scanPosition);
                
            //if(scanDistance > maxSqrDist)
            //    continue;
                
            bool scanIsIn = isIn(scanPosition / iResolution.xy, noise.z);
            if (scanIsIn != fragIsIn)
            {
                if (scanDistance < squaredDistanceToEdge)
                    squaredDistanceToEdge = scanDistance;
            }
        }
    }
    
    float distanceToEdge = sqrt(squaredDistanceToEdge);

    if (fragIsIn)
    {
        // We add 1.0 here since the edge is on the inside
        // of the boundary.
        distanceToEdge = 1.0-distanceToEdge;
    }
        
    distanceToEdge /= halfRange * 2.;
        
    distanceToEdge = 0.5 - distanceToEdge;
    
    distanceToEdge = smoothstep(0., 1., distanceToEdge);
    //distanceToEdge = smoothstep(0., 1., distanceToEdge);
    
    fragColor = vec4(distanceToEdge, distanceToEdge, distanceToEdge, 1.0);

    vec4 oldColor = texture(iChannel2, uv);
    
    // This block of code was to make the terrain
    // "grow" instead of expose the poorly resolved
    // initial state with few samples.
    if(oldColor.a == 0.)
    {
        oldColor.rgb = vec3(25.);
        oldColor.a = 50.;
    }
    fragColor += oldColor;
}

```

## Buffer B

```glsl
float perlin(vec2 uv)
{
    uv /= 4.0;
    vec2 occ = vec2(0);
    float a = 1.0;
    for(int i = 0; i < 6; i++)
    {
        occ += vec2(texture(iChannel0, uv).r, 1) * a;
        uv *= 0.5;
        a *= 2.0;
    }

    float v = occ.x / occ.y;
    
    v = v * 2.0 - 1.0;
    
    v = tanh(v * 2.0);
    
    v = v * 0.5 + 0.5;
    
    //v *= v;
    
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
float maxHeight()
{
    return 200.0 / iChannelResolution[0].x;        
}

float precis()
{
    return 1.5 / iChannelResolution[0].x;        
}


vec3 worldToTerrain(vec3 p)
{
    p.x /= iChannelResolution[0].x / iChannelResolution[0].y;
    p.xz += 0.5;
    p.y /= maxHeight();
    return p;
}

vec3 terrainToWorld(vec3 p)
{
    p.y *= maxHeight();   
    p.xz -= 0.5;
    p.x *= iChannelResolution[0].x / iChannelResolution[0].y;
    return p;
}

float terrain(vec2 p){

    p.x /= iChannelResolution[0].x / iChannelResolution[0].y;
    
    p += 0.5;
    

    
    if(clamp(p,0.0,1.0) != p) return 0.5;
    
    vec4 t = texture(iChannel0, p);
    return t.r / t.a * maxHeight();
}

vec2 calculate_curvature(vec3 p)
{    
    float height_mod = 1.;
    float prec = precis();
    float heightC = terrain(p.xz)*height_mod;
    
    vec3 posC = vec3(p.x, heightC, p.z);
    vec3 posR = p + vec3(prec, 0.0,  0.0);
    vec3 posL = p - vec3(prec, 0.0,  0.0);
    vec3 posT = p + vec3( 0.0, 0.0, prec);
    vec3 posB = p - vec3( 0.0, 0.0, prec);
    
    posR.y = terrain(posR.xz)*height_mod;
    posL.y = terrain(posL.xz)*height_mod;
    posT.y = terrain(posT.xz)*height_mod;
    posB.y = terrain(posB.xz)*height_mod;
    
    vec3 dx = posR - posL;
    vec3 dy = posT - posB;
    
    vec3 normal = normalize(cross(dx, dy));
    
    float curveX = -dot(posC + posC - posR - posL, normal);
    float curveY = -dot(posC + posC - posT - posB, normal);
    
    return vec2(curveX, curveY) / prec;
}

vec3 calculate_normal(vec3 p)
{    
    float height_mod = 1.;
    float prec = precis();
    float heightC = terrain(p.xz)*height_mod;
    
    vec3 posC = vec3(p.x, heightC, p.z);
    vec3 posR = p + vec3(prec, 0.0,  0.0);
    vec3 posL = p - vec3(prec, 0.0,  0.0);
    vec3 posT = p + vec3( 0.0, 0.0, prec);
    vec3 posB = p - vec3( 0.0, 0.0, prec);
    
    posR.y = terrain(posR.xz)*height_mod;
    posL.y = terrain(posL.xz)*height_mod;
    posT.y = terrain(posT.xz)*height_mod;
    posB.y = terrain(posB.xz)*height_mod;
    
    vec3 dx = posR - posL;
    vec3 dy = posT - posB;
    
    return -normalize(cross(dx, dy));
}

bool pointTerrain(vec3 p){
    return terrain(p.xz) >= p.y;
}

vec3 skybox(vec3 dir)
{
    float gradient = dir.y * 0.5 + 0.5;
    
    gradient *= gradient;
    gradient *= gradient;
   
    gradient = 1.-gradient;
    gradient *= gradient;
    gradient *= gradient;
    gradient *= gradient;
    gradient *= gradient;
    gradient = 1.-gradient;
    
    gradient = smoothstep(0., 1., gradient);
    gradient = smoothstep(0., 1., gradient);
    gradient = smoothstep(0., 1., gradient);
    gradient = smoothstep(0., 1., gradient);
    
    vec3 gradient3 = pow(vec3(gradient), vec3(8.0, 1.0, 1.0));
    gradient3 = vec3(1.)-gradient3;
    gradient3 = pow(gradient3, vec3(0.4, 0.5, 4.0));
    gradient3 = vec3(1.)-gradient3;
    
    vec3 color = mix(vec3(0.99,0.99,0.99), vec3(0.01,0.02,0.2), gradient3) * 2.0;
    
    
    return color;
}

vec3 skyboxBlurry(vec3 dir)
{
    float gradient = dir.y * 0.5 + 0.5;
    
    
    
   
    //gradient *= gradient;
   ;
    gradient = 1.-gradient;
    gradient *= gradient * gradient;
    gradient = 1.-gradient;
    
    
    vec3 gradient3 = pow(vec3(gradient), vec3(8.0, 1.0, 1.0));
    gradient3 = vec3(1.)-gradient3;
    gradient3 = pow(gradient3, vec3(0.4, 0.5, 4.0));
    gradient3 = vec3(1.)-gradient3;
    
    vec3 color = mix(vec3(0.99,0.99,0.99), vec3(0.01,0.02,0.2), gradient3) * 2.0;
    
    
    return color;
}

vec3 castRayTerrain2(vec3 camPos, vec3 camDir){
    bounds b;
    b.mini = terrainToWorld(vec3(0));
    b.maxi = terrainToWorld(vec3(1));
    
    vec3 near, far;
    if(!rayIntersectBounds(b, camPos, camDir, near, far))
    {
        return skybox(camDir);
    }
    
    vec3 p = near;
    float minD = 0.0, maxD = 1.0;
    bool hit = false;
    for(float i=0.; i<=1.; i+=.005)
    {
        float j = i*i;
    	p = mix(near, far, j);
        
        maxD = j;
        
        
        if(pointTerrain(p))
        {
            hit = true;
            break;
        }
        
        minD = j;
    }
    
    if(!hit) return skybox(camDir);

    maxD = (minD + maxD) * 0.5;
    float stepSize = (maxD - minD) * 0.5;

    for(int i = 0; i < 5; i++, stepSize *= 0.5)
    {
    	p = mix(near, far, maxD);
        if(!pointTerrain(p))
        {
            minD = maxD;
            maxD += stepSize;
        }
        else
        {
            maxD -= stepSize;
        }
    }
    
    
    vec3 normal = calculate_normal(p);
    vec2 curve2 = calculate_curvature(p);
    
    vec2 posCurve2 = max(vec2(0), curve2);
    vec2 negCurve2 = -min(vec2(0), curve2);
    
    float posCurve = posCurve2.x + posCurve2.y;
    posCurve = posCurve / (0.2 + posCurve);
    float negCurve = negCurve2.x + negCurve2.y;
    negCurve = negCurve / (0.2 + negCurve);
    float curve = (curve2.x + curve2.y);
    curve = curve / (0.2 + abs(curve));
    
    p = vec3(p.x, terrain(p.xz), p.z);

    vec3 lightDir = normalize(vec3(0, 3, 5));
    vec3 lightColor = max(0.0,dot(normal, lightDir)) * vec3(1.0, 0.75, 0.5) * 2.0;
    
    vec3 ambientDir = vec3(0, 1, 0);
    float ambient = dot(normal, ambientDir) * 0.5 + 0.5;
    ambient *= ambient;
    //ambient *= ambient;
    //ambient *= ambient;
    ambient *= curve * 0.5 + 0.5;
    vec3 ambientColor = ambient * skyboxBlurry(normal);
    
    vec3 totalDiffuse = ambientColor + lightColor;

    vec3 eyeDir = normalize(camPos-p);
    vec3 reflec = reflect(-lightDir, normal);
    float spec = max(0.0,dot(eyeDir, reflec));
    spec = pow(spec*1.01, 100.);// * 10.;

    vec3 albedo  = vec3(1.);
    
    float slopeFactor;
    float heightFactor;
    {
        heightFactor = p.y / maxHeight();
        
        heightFactor = smoothstep(0., 1., heightFactor);        
        heightFactor = smoothstep(0., 1., heightFactor);
        
        slopeFactor = abs(normal.y);
        slopeFactor *= slopeFactor;
        slopeFactor *= slopeFactor;
        slopeFactor *= slopeFactor;
        //slopeFactor *= slopeFactor;
        slopeFactor = 1.-slopeFactor;
        slopeFactor *= slopeFactor;
        slopeFactor *= slopeFactor;
        //slopeFactor *= slopeFactor;
        slopeFactor = 1.-slopeFactor;
        //slopeFactor = smoothstep(0., 1., slopeFactor);
        //slopeFactor = smoothstep(0., 1., slopeFactor);
        //slopeFactor = smoothstep(0., 1., slopeFactor);
        
        vec3 flatLow = vec3(0.1, 0.4, 0.05);
        vec3 flatHigh = vec3(0.7, 0.7, 0.5);
        
        vec3 slopeLow = vec3(0.4, 0.7, 0.2);
        vec3 slopeHigh = vec3(0.3, 0.2, 0.05);
        
        vec3 flatColor = mix(flatLow, flatHigh, heightFactor);
        vec3 slopeColor = mix(slopeLow, slopeHigh, heightFactor);
        
        albedo = mix(slopeColor, flatColor, slopeFactor);
    }
    
    albedo = mix(albedo, (albedo*1.0+vec3(0.75))*vec3(0.7,0.65,0.3), posCurve); 
    
    // Add rivers.
    float rivers;
    {
        rivers = max(0.,negCurve - posCurve);
        
        rivers = pow(rivers, pow(2., mix(0.5, -0.7, slopeFactor) + mix(-0.7, 0.7, heightFactor)));   

        //rivers = pow(rivers, heightFactor + 1.0);   
        //rivers = pow(rivers, heightFactor*3.+0.5);         
        
        rivers = smoothstep(0., 1., rivers);
        rivers = smoothstep(0., 1., rivers);
        
        //rivers = 1.-rivers;
        //rivers *= rivers;
        //rivers = 1.-rivers;
        
        albedo = mix(albedo, vec3(0.1,0.2,0.5), rivers); 
    }

    albedo *= albedo;

    //albedo = mix(albedo, curve > 0.0 ? vec3(1,0.7,0.3) : vec3(0.3,0.7,1), vec3(abs(curve))); 

    float specMult = 0.0;
    specMult = mix(specMult, 0.5, posCurve);//length(color)*length(color)*0.2;
    specMult = mix(specMult, 3., rivers);//length(color)*length(color)*0.2;

    vec3 color = albedo*totalDiffuse + lightColor*spec*specMult;
    
    float fogDist = dot(p - camPos, p - camPos);
    
    color = mix(skybox(camDir), color, pow(vec3(0.9), fogDist * vec3(0.2,0.8,4.))); 
    
    return color;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 camRay = normalize(vec3((fragCoord - iResolution.xy * 0.5) / iResolution.yy, 1));
    vec3 camPos = vec3(0.0, maxHeight() * 1.0, -0.5);
    
    vec3 lookDir = vec3(0, maxHeight() * 0.25, 0) - camPos;
    
    quaternion pan = FromAngleAxis(vec3(0, iTime * 0.25, 0));
    quaternion tilt = FromToRotation(vec3(0,0,1), lookDir);
    
    camPos = mul(pan, camPos);
    //camRay = mul(, camRay);
    camRay = mul(mul(pan, tilt), camRay);
    
    vec3 col = castRayTerrain2(camPos, camRay);
    
    col = pow(col, vec3(1./2.2));
    
    col = col / sqrt(0.5+col*col);

    fragColor = vec4(col,1.0);
}
```
