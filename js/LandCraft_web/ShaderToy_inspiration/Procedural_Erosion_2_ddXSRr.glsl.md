
## Common

```glsl
///////////////////////////////////////////////
// HELPER FUNCTIONS
///////////////////////////////////////////////

float dot2(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }
float dot2(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
float dot2(vec4 a, vec4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
int dot2(ivec2 a, ivec2 b) { return a.x * b.x + a.y * b.y; }
int dot2(ivec3 a, ivec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
int dot2(ivec4 a, ivec4 b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }

#define decl_lengthSqr(retType, type) retType lengthSqr(type a) { type t = a; return dot2(t,t); }
decl_lengthSqr(float, vec2)
decl_lengthSqr(float, vec3)
decl_lengthSqr(float, vec4)
decl_lengthSqr(int, ivec2)
decl_lengthSqr(int, ivec3)
decl_lengthSqr(int, ivec4)

#define decl_distSqr(retType, type) retType distSqr(type a, type b) { type diff = a-b; return lengthSqr(diff); }
decl_distSqr(float, vec2)
decl_distSqr(float, vec3)
decl_distSqr(float, vec4)
decl_distSqr(int, ivec2)
decl_distSqr(int, ivec3)
decl_distSqr(int, ivec4)

// I don't trust int(floor(x)).
// So, for positive values:
//
// x = floor(x)
//
// Add 0.5 in case floor(x) rounded it to some floating point garbage like 23.99999
// x += 0.5
//
// ix = int(x) // THEN cast to an int.
//
// Similar thing happens for negative values.
#define decl_floorToInt(retType, type) retType floorToInt(type a) { return retType(floor(a) + sign(a) * type(0.5)); }
decl_floorToInt(int, float)
decl_floorToInt(ivec2, vec2)
decl_floorToInt(ivec3, vec3)
decl_floorToInt(ivec4, vec4)
   
   // RNG
uint wang_hash(inout uint seed)
{
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}

float RandomFloat01(inout uint state)
{
    return float(wang_hash(state)) / 4294967296.0;
}

// Credit: https://www.shadertoy.com/view/DsBGzy by sh1boot
vec3 hsv2rgb(vec3 hsv) 
{
    vec3 h3 = mod(6.0 * hsv.x + vec3(5.0, 3.0, 1.0), 6.0);
    h3 = min(h3, 4.0 - h3);
    h3 = clamp(h3, 0.0, 1.0);
    return hsv.z - hsv.z * hsv.y * h3;
}

// Credit: https://www.shadertoy.com/view/DsBGzy by sh1boot
vec3 palette(int i) 
{
    float f = float(i);
    float h =  mod(1.618033988749894848204586834 * f, 1.0);
    float s = exp(-0.00025 * f) * 0.65 + 0.25;
    float v = 1.0;
    return hsv2rgb(vec3(h, s, v));
}

vec4 cube(samplerCube cube, vec2 uv)
{
    uv -= floor(uv);
    uv = uv * 2.0 - 1.0;
    
    vec3 ray = vec3(uv, 1);
    
    return texture(cube, ray);
}

vec4 cubeFetch(samplerCube cube, ivec2 coord, int mip)
{
    vec2 uv = (vec2(coord) + vec2(0.5)) / vec2(textureSize(cube, mip));
    //uv = clamp(uv, vec2(0), vec2(1));
    uv -= floor(uv);
    uv = uv * 2.0 - 1.0;
    
    vec3 ray = vec3(uv, 1);
    
    return textureLod(cube, ray, float(mip));
}

vec4 cubeLod(samplerCube cube, vec2 uv, float mip)
{
    //uv = clamp(uv, vec2(0), vec2(1));
    uv -= floor(uv);
    uv = uv * 2.0 - 1.0;
    
    vec3 ray = vec3(uv, 1);
    
    return textureLod(cube, ray, mip);
}

// Drawing helpers --------------------------------

void blend(in vec4 src, inout vec4 dest)
{
    dest.rgb = mix(dest.rgb, src.rgb, src.a);
    dest.a = mix(dest.a, 1.0, src.a);    
}

void drawCircle(vec2 center, vec2 frag, float radius, float lineWidth, vec4 lineColor, inout vec4 color)
{
    float distanceToEdge = abs(radius - distance(frag, center));

    float circle = smoothstep
    (
        0.0, 
        1.0,
        1.0 - (distanceToEdge / lineWidth) 
    );
                
    lineColor.a *= circle;
    blend(lineColor, color);
}

void drawSDF(vec3 diff_dist, inout vec4 color)
{
    float d = diff_dist.z;
    vec2 grad = diff_dist.xy/(d);
    
    vec3 c = normalize(vec3(grad,sign(d))) * 0.5 + 0.5;
	c *= 1. - exp2(-12. * abs(d));
	c *= .8 + .2 * cos(120.*d);

    color.rgb = c;
}

vec3 unsignedValueToColor(float value)
{
    float valueLog = log2(value + 0.0625);
    valueLog /= 16.0;
    //value /= 1. + abs(valueLog);
    valueLog = valueLog * 0.5 + 0.5;
    
    value *= 10.0;
    value += 0.125;
    value /= 1. + abs(value);
    
    return hsv2rgb(vec3(valueLog, 1, value));
}

// Font helpers ------------------------------
#define FONTSAMPLERSIZE vec2(1024, 1024)
#define FONTSAMPLERSIZEI ivec2(1024, 1024)
#define FONTELEMENTCOUNT  vec2(16, 16)
#define FONTELEMENTCOUNTI  ivec2(16, 16)
#define FONTELEMENTSIZE  vec2(64, 64)
#define FONTELEMENTSIZEI  ivec2(64, 64)
vec4 sampleFontElement(sampler2D fontSampler, in vec2 fragCoord, in ivec2 element)
{
    vec2 samplerCoordMin = vec2(element * FONTELEMENTSIZEI); 
    vec2 samplerCoordMax = vec2(element * FONTELEMENTSIZEI + FONTELEMENTSIZEI); 
    
    vec2 elementUV = (fragCoord - samplerCoordMin) / FONTELEMENTSIZE;
    elementUV -= floor(elementUV);
    
    vec2 samplerCoord = mix(samplerCoordMin, samplerCoordMax, elementUV);
    vec2 samplerUV = samplerCoord / FONTSAMPLERSIZE;
    
    return texture(fontSampler, samplerUV);
}

vec4 sampleFontElementColor(sampler2D fontSampler, in vec2 fragCoord, in ivec2 element, in vec4 color)
{
    float opacity = sampleFontElement(fontSampler, fragCoord, element).r;
    
    color.a *= opacity;
    
    return color;
}

vec3 getNoise(int iFrame)
{
    vec3 noise = vec3(0);

    // Calculate oising over time.
    vec3 temporalNoise;
    {
        // Use the golden ratio as it should land
        // on all fractional values eventually.
        temporalNoise = vec3(iFrame, iFrame+1, iFrame+2);
        temporalNoise *= vec3(0.7548776662, 0.56984029, 0.618033988749);

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
    
    return noise;
}

/////////////////////////////////////////////////////////////
// C# Code used to precompute coordinate lookup tables:
/*
using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
					
public class Program
{
	static void GetSubsampleCoords(int index, out int[] c0, out int[] c1, out int[] c2, out int[] c3)
	{
		int r = 0x1B; // 00 10 01 11;

		if((index & 1) != 0)
		{
			const int xor = 0xAA; // 10 10 10 10
			r = r ^ xor;
		}

		if((index & 2) != 0)
		{
			const int xor = 0x55; // 01 01 01 01
			r = r ^ xor;
		}

		if((index & 4) != 0)
		{
			const int xor = 0x3C; // 00 11 11 00
			r = r ^ xor;
		}

		int[] c0_m = new int[2]
		{
			0x80, // 10 00 00 00
			0x40 // 01 00 00 00
		};

		int[] c0_s = new int[2]
		{
			7,
			6
		};
		
		int[] c1_m = new int[2]
		{
			0x20, // 00 10 00 00
			0x10  // 00 01 00 00
		};

		int[] c1_s = new int[2]
		{
			5, 
			4 
		};

		int[] c2_m = new int[2]
		{
			0x08, // 00 00 10 00
			0x04  // 00 00 01 00
		};

		int[] c2_s = new int[2]
		{
			3, 
			2
		};
		
		int[] c3_m = new int[2]
		{
			0x02, // 00 00 00 10
			0x01  // 00 00 00 01
		};

		int[] c3_s = new int[2]
		{
			1,
			0
		};
		
		c0 = new int[2];		
		c1 = new int[2];
		c2 = new int[2];
		c3 = new int[2];
		for(int i = 0; i < 2; i++)
		{
			c0[i] = (r & c0_m[i]) >> c0_s[i];
			c1[i] = (r & c1_m[i]) >> c1_s[i];
			c2[i] = (r & c2_m[i]) >> c2_s[i];
			c3[i] = (r & c3_m[i]) >> c3_s[i];
		}
	}
	
	static string AsString(int[] v)
	{
		string s = "ivec2(";
		
		for(int i = 0; i < v.Length; i++)
		{
			if(i > 0)
			{
				s = s + ", ";
			}
			s = s + v[i].ToString();
		}
		
		s += "), ";
		return s;
	}
	
	static Vector2 GetTestCoord(int index)
	{
		Vector2 v = default(Vector2);
		
		// abs(x) > abs(y)?
		if((index & 4) != 0)
		{
			v.X = 0.25f;
			v.Y = 0.125f;
		}	
		else
		{
			v.X = 0.125f;
			v.Y = 0.25f;
		}
		
		// x < 0?
		if((index & 1) != 0)
		{
			v.X = -v.X;
		}
		
		// y < 0?
		if((index & 2) != 0)
		{
			v.Y = -v.Y;
		}
		
		return -v;
	}
	
	static string AsString(Vector2 v)
	{
		return string.Format
		(
			"ivec2({0}{1}, {2}{3})", 
			 v.X < 0 ? '-' : (v.X > 0 ? '+' : ' '), 
			 Math.Abs(v.X), 
			 v.Y < 0 ? '-' : (v.Y > 0 ? '+' : ' '), 
			 Math.Abs(v.Y)
	     );
	}
	
	static string AsString(Vector2[] v)
	{
		string s = "";
		
		for(int i = 0; i < v.Length; i++)
		{
			if(i > 0)
			{
				s = s + ", ";
				
				if(i % 4 == 0)
					s = s + '\n';
			}
			s = s + AsString(v[i]);//.ToString();
		}
		
		s += ",";
		return s;
	}
	
	public static void Main()
	{
		for(int i = 0; i < 8; i++)
		{
			int[] a, b, c, d;
			GetSubsampleCoords(i, out a, out b, out c, out d);
			
			
			Console.WriteLine(AsString(a) 
							+ AsString(b)
							+ AsString(c)
							+ AsString(d));
		}
		
		Console.WriteLine("________________");
		Console.WriteLine("                ");
		
		const int radius = 4;
		const int diameter = radius+radius+1;
		const int capacity = diameter * diameter - 1;
		List<Vector2> pixelCoords = new List<Vector2>(capacity);
		
		const float maxDistanceSqr = (radius + 0.5f) * (radius + 0.5f);

		for(int x = -radius; x <= radius; x++)
		{
			for(int y = -radius; y <= radius; y++)
			{
				if(x == 0 && y == 0) continue;
				pixelCoords.Add(new Vector2(x,y));
			}
		}
		
		for(int i = 0; i < 8; i++)
		{
			Vector2 testCoord = GetTestCoord(i);
			
			var pixelsSorted = pixelCoords.Where(p => Vector2.DistanceSquared(p, testCoord) <= maxDistanceSqr).OrderBy(p => Vector2.DistanceSquared(p, testCoord)).ToArray();
			
			//Console.WriteLine(pixelsSorted.Length);
			Console.WriteLine(AsString(pixelsSorted));
			Console.WriteLine("");
		}
	}
}
*/
```

## Buffer A

```glsl
// HOLY ---- I DID IT
// Broad crawl outward, then a recursive mipmap quadtree search.
// FAST runtime, FAST compilation (comparitavely).
// Should compile in about 7.2 seconds, 144+ FPS @ 900x506

// Controls =================================
// Frag + Click: Change the distance checked position.

// Colors ===================================
// Black Cell:       Cell fully empty.
// Blue Cell:        Cell partially occupied (contains an edge).
// Light Blue Cell:  Cell fully occupied.
// Red-Tinted Cells: Cells checked for occupancy.
// Circle:           (Color) Gradient (Radius) Final distance.

const int mipCount = 8;
const int maxMip = mipCount - 1;

bool containsEdge(bool mySign, ivec2 coord, int mipLevel)
{
    float occupancy = cubeFetch(iChannel0, coord, mipLevel).r;

    if(mySign)
    {
        return occupancy < 1.;
    }
    else
    {
        return occupancy > 0.;
    }
}

bool isOccupied(vec2 uv)
{
    float occupancy = cubeLod(iChannel0, uv, 0.0).r;
    return occupancy > 0.;
}

const int SampleCoordCount = 64*8;
const ivec2 SampleCoords[64*8] = ivec2[64*8]
(
    ivec2( 0, -1), ivec2(-1,  0), ivec2(-1, -1), ivec2(+1,  0), 
    ivec2( 0, +1), ivec2(+1, -1), ivec2(-1, +1), ivec2(+1, +1), 
    ivec2( 0, -2), ivec2(-2,  0), ivec2(-1, -2), ivec2(-2, -1), 
    ivec2(+1, -2), ivec2(+2,  0), ivec2(-2, +1), ivec2( 0, +2), 
    ivec2(+2, -1), ivec2(-1, +2), ivec2(+2, +1), ivec2(+1, +2), 
    ivec2(-2, -2), ivec2( 0, -3), ivec2(+2, -2), ivec2(-3,  0), 
    ivec2(-1, -3), ivec2(-2, +2), ivec2(-3, -1), ivec2(+1, -3), 
    ivec2(+2, +2), ivec2(-3, +1), ivec2(+3,  0), ivec2(+3, -1), 
    ivec2( 0, +3), ivec2(-2, -3), ivec2(-3, -2), ivec2(-1, +3), 
    ivec2(+3, +1), ivec2(+1, +3), ivec2(+2, -3), ivec2(+3, -2), 
    ivec2(-3, +2), ivec2(-2, +3), ivec2( 0, -4), ivec2(-1, -4), 
    ivec2(+3, +2), ivec2(-4,  0), ivec2(+2, +3), ivec2(+1, -4), 
    ivec2(-4, -1), ivec2(-3, -3), ivec2(-4, +1), ivec2(+4,  0), 
    ivec2(+3, -3), ivec2(-2, -4), ivec2(+4, -1), ivec2(-4, -2), 
    ivec2( 0, +4), ivec2(+2, -4), ivec2(+4, +1), ivec2(-3, +3), 
    ivec2(-1, +4), ivec2(+1, +4), ivec2(-4, +2), ivec2(+4, -2),

    ivec2( 0, -1), ivec2(+1,  0), ivec2(-1,  0), ivec2(+1, -1), 
    ivec2( 0, +1), ivec2(-1, -1), ivec2(+1, +1), ivec2(-1, +1), 
    ivec2( 0, -2), ivec2(+2,  0), ivec2(+1, -2), ivec2(+2, -1), 
    ivec2(-1, -2), ivec2(-2,  0), ivec2(-2, -1), ivec2( 0, +2), 
    ivec2(+2, +1), ivec2(+1, +2), ivec2(-2, +1), ivec2(-1, +2), 
    ivec2(+2, -2), ivec2(-2, -2), ivec2( 0, -3), ivec2(+1, -3), 
    ivec2(+3,  0), ivec2(+2, +2), ivec2(-1, -3), ivec2(+3, -1), 
    ivec2(-2, +2), ivec2(-3,  0), ivec2(+3, +1), ivec2(-3, -1), 
    ivec2( 0, +3), ivec2(+2, -3), ivec2(-3, +1), ivec2(+1, +3), 
    ivec2(+3, -2), ivec2(-1, +3), ivec2(-2, -3), ivec2(-3, -2), 
    ivec2(+3, +2), ivec2( 0, -4), ivec2(+2, +3), ivec2(-3, +2), 
    ivec2(+1, -4), ivec2(-2, +3), ivec2(+4,  0), ivec2(-1, -4), 
    ivec2(+4, -1), ivec2(+3, -3), ivec2(+4, +1), ivec2(-4,  0), 
    ivec2(-3, -3), ivec2(-4, -1), ivec2(+2, -4), ivec2( 0, +4), 
    ivec2(+4, -2), ivec2(-4, +1), ivec2(-2, -4), ivec2(+1, +4), 
    ivec2(+3, +3), ivec2(-1, +4), ivec2(-4, -2), ivec2(+4, +2),

    ivec2( 0, +1), ivec2(-1,  0), ivec2(-1, +1), ivec2(+1,  0), 
    ivec2( 0, -1), ivec2(+1, +1), ivec2(-1, -1), ivec2(+1, -1), 
    ivec2( 0, +2), ivec2(-2,  0), ivec2(-1, +2), ivec2(-2, +1), 
    ivec2(+1, +2), ivec2(+2,  0), ivec2(-2, -1), ivec2( 0, -2), 
    ivec2(+2, +1), ivec2(-1, -2), ivec2(+2, -1), ivec2(+1, -2), 
    ivec2(-2, +2), ivec2( 0, +3), ivec2(+2, +2), ivec2(-3,  0), 
    ivec2(-1, +3), ivec2(-2, -2), ivec2(-3, +1), ivec2(+1, +3), 
    ivec2(+2, -2), ivec2(-3, -1), ivec2(+3,  0), ivec2(+3, +1), 
    ivec2( 0, -3), ivec2(-2, +3), ivec2(-3, +2), ivec2(-1, -3), 
    ivec2(+3, -1), ivec2(+1, -3), ivec2(+2, +3), ivec2(+3, +2), 
    ivec2(-3, -2), ivec2(-2, -3), ivec2( 0, +4), ivec2(-1, +4), 
    ivec2(+3, -2), ivec2(-4,  0), ivec2(+2, -3), ivec2(+1, +4), 
    ivec2(-4, +1), ivec2(-3, +3), ivec2(-4, -1), ivec2(+4,  0), 
    ivec2(+3, +3), ivec2(-2, +4), ivec2(+4, +1), ivec2(-4, +2), 
    ivec2( 0, -4), ivec2(+2, +4), ivec2(+4, -1), ivec2(-3, -3), 
    ivec2(-1, -4), ivec2(+1, -4), ivec2(-4, -2), ivec2(+4, +2),

    ivec2( 0, +1), ivec2(+1,  0), ivec2(-1,  0), ivec2(+1, +1), 
    ivec2( 0, -1), ivec2(-1, +1), ivec2(+1, -1), ivec2(-1, -1), 
    ivec2( 0, +2), ivec2(+2,  0), ivec2(+1, +2), ivec2(+2, +1), 
    ivec2(-1, +2), ivec2(-2,  0), ivec2(-2, +1), ivec2( 0, -2), 
    ivec2(+2, -1), ivec2(+1, -2), ivec2(-2, -1), ivec2(-1, -2), 
    ivec2(+2, +2), ivec2(-2, +2), ivec2( 0, +3), ivec2(+1, +3), 
    ivec2(+3,  0), ivec2(+2, -2), ivec2(-1, +3), ivec2(+3, +1), 
    ivec2(-2, -2), ivec2(-3,  0), ivec2(+3, -1), ivec2(-3, +1), 
    ivec2( 0, -3), ivec2(+2, +3), ivec2(-3, -1), ivec2(+1, -3), 
    ivec2(+3, +2), ivec2(-1, -3), ivec2(-2, +3), ivec2(-3, +2), 
    ivec2(+3, -2), ivec2( 0, +4), ivec2(+2, -3), ivec2(-3, -2), 
    ivec2(+1, +4), ivec2(-2, -3), ivec2(+4,  0), ivec2(-1, +4), 
    ivec2(+4, +1), ivec2(+3, +3), ivec2(+4, -1), ivec2(-4,  0), 
    ivec2(-3, +3), ivec2(-4, +1), ivec2(+2, +4), ivec2( 0, -4), 
    ivec2(+4, +2), ivec2(-4, -1), ivec2(-2, +4), ivec2(+1, -4), 
    ivec2(+3, -3), ivec2(-1, -4), ivec2(-4, +2), ivec2(+4, -2),

    ivec2(-1,  0), ivec2( 0, -1), ivec2(-1, -1), ivec2( 0, +1), 
    ivec2(+1,  0), ivec2(-1, +1), ivec2(+1, -1), ivec2(+1, +1), 
    ivec2(-2,  0), ivec2( 0, -2), ivec2(-2, -1), ivec2(-1, -2), 
    ivec2(-2, +1), ivec2( 0, +2), ivec2(-1, +2), ivec2(+1, -2), 
    ivec2(+2,  0), ivec2(+2, -1), ivec2(+1, +2), ivec2(+2, +1), 
    ivec2(-2, -2), ivec2(-3,  0), ivec2(-2, +2), ivec2(-3, -1), 
    ivec2( 0, -3), ivec2(+2, -2), ivec2(-3, +1), ivec2(-1, -3), 
    ivec2(+2, +2), ivec2( 0, +3), ivec2(+1, -3), ivec2(-1, +3), 
    ivec2(+3,  0), ivec2(-3, -2), ivec2(-2, -3), ivec2(+1, +3), 
    ivec2(+3, -1), ivec2(+3, +1), ivec2(-3, +2), ivec2(-2, +3), 
    ivec2(+2, -3), ivec2(-4,  0), ivec2(+3, -2), ivec2(-4, -1), 
    ivec2(+2, +3), ivec2( 0, -4), ivec2(+3, +2), ivec2(-4, +1), 
    ivec2(-1, -4), ivec2(-3, -3), ivec2(+1, -4), ivec2( 0, +4), 
    ivec2(-3, +3), ivec2(-4, -2), ivec2(-1, +4), ivec2(-2, -4), 
    ivec2(+4,  0), ivec2(-4, +2), ivec2(+1, +4), ivec2(+3, -3), 
    ivec2(+4, -1), ivec2(+4, +1), ivec2(-2, +4), ivec2(+2, -4),

    ivec2(+1,  0), ivec2( 0, -1), ivec2( 0, +1), ivec2(+1, -1), 
    ivec2(-1,  0), ivec2(+1, +1), ivec2(-1, -1), ivec2(-1, +1), 
    ivec2(+2,  0), ivec2( 0, -2), ivec2(+2, -1), ivec2(+1, -2), 
    ivec2(+2, +1), ivec2( 0, +2), ivec2(-2,  0), ivec2(-1, -2), 
    ivec2(+1, +2), ivec2(-2, -1), ivec2(-1, +2), ivec2(-2, +1), 
    ivec2(+2, -2), ivec2(+2, +2), ivec2(+3,  0), ivec2( 0, -3), 
    ivec2(+3, -1), ivec2(-2, -2), ivec2(+1, -3), ivec2(+3, +1), 
    ivec2(-2, +2), ivec2(-1, -3), ivec2( 0, +3), ivec2(+1, +3), 
    ivec2(-3,  0), ivec2(+3, -2), ivec2(-3, -1), ivec2(-1, +3), 
    ivec2(+2, -3), ivec2(-3, +1), ivec2(+3, +2), ivec2(+2, +3), 
    ivec2(-2, -3), ivec2(-3, -2), ivec2(+4,  0), ivec2(-2, +3), 
    ivec2(+4, -1), ivec2(-3, +2), ivec2( 0, -4), ivec2(+4, +1), 
    ivec2(+1, -4), ivec2(+3, -3), ivec2(-1, -4), ivec2( 0, +4), 
    ivec2(+3, +3), ivec2(+1, +4), ivec2(+4, -2), ivec2(-4,  0), 
    ivec2(+2, -4), ivec2(-1, +4), ivec2(+4, +2), ivec2(-4, -1), 
    ivec2(-3, -3), ivec2(-4, +1), ivec2(-2, -4), ivec2(+2, +4),

    ivec2(-1,  0), ivec2( 0, +1), ivec2(-1, +1), ivec2( 0, -1), 
    ivec2(+1,  0), ivec2(-1, -1), ivec2(+1, +1), ivec2(+1, -1), 
    ivec2(-2,  0), ivec2( 0, +2), ivec2(-2, +1), ivec2(-1, +2), 
    ivec2(-2, -1), ivec2( 0, -2), ivec2(-1, -2), ivec2(+1, +2), 
    ivec2(+2,  0), ivec2(+2, +1), ivec2(+1, -2), ivec2(+2, -1), 
    ivec2(-2, +2), ivec2(-3,  0), ivec2(-2, -2), ivec2(-3, +1), 
    ivec2( 0, +3), ivec2(+2, +2), ivec2(-3, -1), ivec2(-1, +3), 
    ivec2(+2, -2), ivec2( 0, -3), ivec2(+1, +3), ivec2(-1, -3), 
    ivec2(+3,  0), ivec2(-3, +2), ivec2(-2, +3), ivec2(+1, -3), 
    ivec2(+3, +1), ivec2(+3, -1), ivec2(-3, -2), ivec2(-2, -3), 
    ivec2(+2, +3), ivec2(-4,  0), ivec2(+3, +2), ivec2(-4, +1), 
    ivec2(+2, -3), ivec2( 0, +4), ivec2(+3, -2), ivec2(-4, -1), 
    ivec2(-1, +4), ivec2(-3, +3), ivec2(+1, +4), ivec2( 0, -4), 
    ivec2(-3, -3), ivec2(-4, +2), ivec2(-1, -4), ivec2(-2, +4), 
    ivec2(+4,  0), ivec2(-4, -2), ivec2(+1, -4), ivec2(+3, +3), 
    ivec2(+4, +1), ivec2(+4, -1), ivec2(-2, -4), ivec2(+2, +4),

    ivec2(+1,  0), ivec2( 0, +1), ivec2( 0, -1), ivec2(+1, +1), 
    ivec2(-1,  0), ivec2(+1, -1), ivec2(-1, +1), ivec2(-1, -1), 
    ivec2(+2,  0), ivec2( 0, +2), ivec2(+2, +1), ivec2(+1, +2), 
    ivec2(+2, -1), ivec2( 0, -2), ivec2(-2,  0), ivec2(-1, +2), 
    ivec2(+1, -2), ivec2(-2, +1), ivec2(-1, -2), ivec2(-2, -1), 
    ivec2(+2, +2), ivec2(+2, -2), ivec2(+3,  0), ivec2( 0, +3), 
    ivec2(+3, +1), ivec2(-2, +2), ivec2(+1, +3), ivec2(+3, -1), 
    ivec2(-2, -2), ivec2(-1, +3), ivec2( 0, -3), ivec2(+1, -3), 
    ivec2(-3,  0), ivec2(+3, +2), ivec2(-3, +1), ivec2(-1, -3), 
    ivec2(+2, +3), ivec2(-3, -1), ivec2(+3, -2), ivec2(+2, -3), 
    ivec2(-2, +3), ivec2(-3, +2), ivec2(+4,  0), ivec2(-2, -3), 
    ivec2(+4, +1), ivec2(-3, -2), ivec2( 0, +4), ivec2(+4, -1), 
    ivec2(+1, +4), ivec2(+3, +3), ivec2(-1, +4), ivec2( 0, -4), 
    ivec2(+3, -3), ivec2(+1, -4), ivec2(+4, +2), ivec2(-4,  0), 
    ivec2(+2, +4), ivec2(-1, -4), ivec2(+4, -2), ivec2(-4, +1), 
    ivec2(-3, +3), ivec2(-4, -1), ivec2(-2, +4), ivec2(+2, -4)
);

const ivec2 SubsampleCoords[4*8] = ivec2[4*8] 
(
    ivec2(0, 0), ivec2(0, 1), ivec2(1, 0), ivec2(1, 1),
    ivec2(1, 0), ivec2(1, 1), ivec2(0, 0), ivec2(0, 1),
    ivec2(0, 1), ivec2(0, 0), ivec2(1, 1), ivec2(1, 0),
    ivec2(1, 1), ivec2(1, 0), ivec2(0, 1), ivec2(0, 0),
    ivec2(0, 0), ivec2(1, 0), ivec2(0, 1), ivec2(1, 1),
    ivec2(1, 0), ivec2(0, 0), ivec2(1, 1), ivec2(0, 1),
    ivec2(0, 1), ivec2(1, 1), ivec2(0, 0), ivec2(1, 0),
    ivec2(1, 1), ivec2(0, 1), ivec2(1, 0), ivec2(0, 0)
);

int GetSampleIndexOffset(vec2 testToSample)
{
    int index;
    index  = testToSample.x < 0. ? (1*64) : 0;
    index |= testToSample.y < 0. ? (2*64) : 0;
    index |= abs(testToSample.x) > abs(testToSample.y) ? (4*64) : 0;
    return index;
}

int GetSubsampleIndexOffset(vec2 testToSample)
{
    int index;
    index  = testToSample.x < 0. ? (1*4) : 0;
    index |= testToSample.y < 0. ? (2*4) : 0;
    index |= abs(testToSample.x) > abs(testToSample.y) ? (4*4) : 0;
    return index;
}

// Returns the minimum distance from
// a given coordinate to a square whose minimum
// is coord, and whose maximum is coord+(1,1).
vec3 get_diff_minDistSqr(vec2 point, ivec2 coord)
{
    vec2 squareCenter = vec2(coord) + vec2(0.5);
    vec2 pointToCenter = squareCenter - point;
    vec2 minOffset = clamp(pointToCenter, vec2(-0.5), vec2(0.5));
    vec2 pointToMin = pointToCenter - minOffset;
    return vec3(pointToMin, lengthSqr(pointToMin));
}

struct SubSampleArgs
{
    ivec3 smpCrd_mip;
    vec3 diff_distSqr;
    vec2 testCoord;
};

#define MIP smpCrd_mip.z
#define SAMPLECOORD smpCrd_mip.xy
#define TESTCOORD testCoord
#define DIFFDISTSQR diff_distSqr

void IncreaseDepth(inout SubSampleArgs a)
{
    a.diff_distSqr *= vec3(2.0,2.0,4.0);
    a.testCoord *= vec2(2.0,2.0);
    a.smpCrd_mip += ivec3(a.smpCrd_mip.xy, -1);      
}

SubSampleArgs GetNextArgs(in SubSampleArgs a, ivec2 coord)
{
    a.SAMPLECOORD = coord;
    a.diff_distSqr *= vec3(2.0,2.0,4.0);
    a.testCoord *= vec2(2.0,2.0);
    a.smpCrd_mip += ivec3(a.smpCrd_mip.xy, -1);    
    return a;
}

// Define the function in a way that recursion is easy to type.
// Too bad multi-line macros are not supported.
#define SSH(name) vec3 name
#define SSA00 (
#define SSA01     SubSampleArgs a,
#define SSA02     in bool    testOccupied
#define SSA03 )
#define SSA04 {    
    
                  // Get the subsample coordinates sorted by distance to the sample.
#define SSA05     int iSubOffset = GetSubsampleIndexOffset(vec2(a.SAMPLECOORD + ivec2(1)) - a.TESTCOORD);   
    
                 // Go through subsamples from closest to furthest.
#define SSA06    for(int iSub = 0; iSub < 4; iSub++)
#define SSA07    {
#define SSA08        ivec2 subsamplingCoord = a.SAMPLECOORD + SubsampleCoords[iSub + iSubOffset];
#define SSA09        vec3 subsampleDiffDistSqr = get_diff_minDistSqr(a.TESTCOORD, subsamplingCoord);
        
#define SSA10        if(!(subsampleDiffDistSqr.z < a.DIFFDISTSQR.z))
#define SSA11            continue;
    
#define SSA12        if(!containsEdge(testOccupied, subsamplingCoord, a.MIP))
#define SSA13            continue;

                     // Recursive   
#define SSR(name)    subsampleDiffDistSqr = name                     
#define SSB00        (
#define SSB01            GetNextArgs(a, subsamplingCoord),
#define SSB02            testOccupied
#define SSB03        );

#define SSC00        if(a.DIFFDISTSQR.z > subsampleDiffDistSqr.z) a.DIFFDISTSQR = subsampleDiffDistSqr;
#define SSC01     }
    
#define SSC02     return a.DIFFDISTSQR * vec3(0.5, 0.5, 0.25);
#define SSC03 }

// Combine all of the lines into one macro.
#define SSA_0 SSA00 SSA01 SSA02 SSA03 SSA04 SSA05 SSA06 SSA07 SSA08 SSA09
#define SSA_1 SSA10 SSA11 SSA12 SSA13 

#define SSB_0 SSB00 SSB01 SSB02 SSB03

#define SSC_0 SSC00 SSC01 SSC02 SSC03

#define SS_HEADER(n)        SSH(n)
#define SS_START            SSA_0 SSA_1
#define SS_RECURSIONCALL(n) SSR(n)
#define SS_RECURSIVE        SSB_0
#define SS_END              SSC_0

#define SS_DECLARE_TERMINATING(name) SS_HEADER(name) SS_START SS_END

#define SS_DECLARE_RECURSIVE(name, calls) SS_HEADER(name) SS_START SS_RECURSIONCALL(calls) SS_RECURSIVE SS_END

#define RECURSION_DEPTH 7
SS_DECLARE_TERMINATING(SubSample7)
SS_DECLARE_RECURSIVE(SubSample6, SubSample7)
SS_DECLARE_RECURSIVE(SubSample5, SubSample6)
SS_DECLARE_RECURSIVE(SubSample4, SubSample5)
SS_DECLARE_RECURSIVE(SubSample3, SubSample4)
SS_DECLARE_RECURSIVE(SubSample2, SubSample3)
SS_DECLARE_RECURSIVE(SubSample, SubSample2)

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 noise = getNoise(iFrame);
    
    ivec2 maxMipSize = textureSize(iChannel0, maxMip); 
    
    // Set up fragment coordinates.
    vec2 fragUV = fragCoord / iChannelResolution[0].xy;
    vec2 fragTestCoord = fragUV * vec2(maxMipSize) 
                       + noise.xy * (1.0 / float(1 << RECURSION_DEPTH));
    ivec2 fragTestIntCoord = floorToInt(fragTestCoord);
    bool fragIsOccupied = isOccupied(fragUV); 
    
    const float maxDist = 2.0;
    vec3 diff_distSqr = vec3(maxDist, 0.0, maxDist*maxDist);

    int iCoordOffset = GetSampleIndexOffset((vec2(fragTestIntCoord) + vec2(0.5) - fragTestCoord));
    for(int iCoord = -1; iCoord < 64; iCoord++)
    {
        ivec2 samplingCoord = fragTestIntCoord;
        if(iCoord >= 0)
        {
            samplingCoord += SampleCoords[iCoord + iCoordOffset];
        }
        
        vec3 sampleDiffDistSqr = get_diff_minDistSqr(fragTestCoord, samplingCoord);
        
        // Need to add a bias for some reason.
        if(sampleDiffDistSqr.z >= diff_distSqr.z + 1.5)
            break;
            
        if(sampleDiffDistSqr.z >= diff_distSqr.z)
            continue; 
            
        if(!containsEdge(fragIsOccupied, samplingCoord, maxMip))
            continue;
            
        SubSampleArgs a;
        a.DIFFDISTSQR = diff_distSqr;
        a.SAMPLECOORD = samplingCoord;
        a.MIP = maxMip;
        a.TESTCOORD = fragTestCoord;
        IncreaseDepth(a);
        diff_distSqr = SubSample
        (
            a,
            fragIsOccupied
        ); 
    }
         
    vec3 diff_dist = vec3(diff_distSqr.xy, sqrt(diff_distSqr.z)) / (fragIsOccupied ? -maxDist : maxDist);     
         
    diff_dist.xy = normalize(diff_dist.xy);
    
    diff_dist.z = diff_dist.z * 0.5 + 0.5;
    diff_dist.z = smoothstep(0., 1., diff_dist.z);
    diff_dist.z = smoothstep(0., 1., diff_dist.z);
    diff_dist.z = smoothstep(0., 1., diff_dist.z);
    //diff_dist.z = smoothstep(0., 1., diff_dist.z);
    
         
    fragColor = vec4(diff_dist, 1.0);
    if(isnan(fragColor.z + fragColor.a)) fragColor = vec4(0);
    if(isinf(fragColor.z + fragColor.a)) fragColor = vec4(0);
    
    // Accumulate samples over time.
    if(iFrame > 1 && iMouse.z <= 0.)
    {
        vec2 uv = (fragCoord.xy) / iChannelResolution[1].xy;
        vec4 oldColor = texture(iChannel1, uv);
        fragColor += oldColor;
    }
    
}

```

## Cube A

```glsl
// This is a really stupid work-around to enforcing power-of-two buffers... but whatever.

float perlin(vec2 uv, inout uint rngState)
{
    uv += vec2(RandomFloat01(rngState), RandomFloat01(rngState));
    uv /= 1.0;
    vec2 occ = vec2(0);
    float a = 1.0;
    for(int i = 0; i <8; i++)
    {
        occ += vec2(texture(iChannel0, uv).r, 1) * a;
        uv *= 0.5;
        a *= 1.7;
        uv += vec2(RandomFloat01(rngState), RandomFloat01(rngState));
    }
    float v = occ.x / occ.y;
    
    return v;
       
}

void mainCubemap( out vec4 fragColor, in vec2 fragCoord, in vec3 rayOri, in vec3 rayDir )
{
    vec2 uv = fragCoord / iResolution.xy; 
    
    vec2 uvTile = uv * 2.0;
    
    if(uvTile.x > 1.0)
        uvTile.x = 2.0-uvTile.x;
    if(uvTile.y > 1.0)
        uvTile.y = 2.0-uvTile.y;
        
    
    uvTile = smoothstep(0., 1., uvTile);
        
    uint rngState = 14114u;
    
    vec2 uv00 = uv + vec2(0.5, 0.5);
    uv00 -= floor(uv00);
    vec2 uv01 = uv + vec2(0.0, 0.5);
    uv01 -= floor(uv01);
    vec2 uv10 = uv + vec2(0.5, 0.0);
    uv10 -= floor(uv10);
    vec2 uv11 = uv + vec2(0.0, 0.0);
    uv11 -= floor(uv11);
    
    float height00 = perlin(uv00, rngState);
    float height01 = perlin(uv01, rngState);
    float height10 = perlin(uv10, rngState);
    float height11 = perlin(uv11, rngState);
    
     //height00 = 0.;
     //height01 = 0.;
     //height10 = 0.;
     //height11 = 0.;
    
    float height0X = mix(height00, height01, uvTile.x);
    float height1X = mix(height10, height11, uvTile.x);
    float v = mix(height0X, height1X, uvTile.y);
    
    //v = height01;
    /*float height1 = perlin(uv, rngState);
    
    vec2 offsetUv = uv;
    offsetUv += vec2(0.5, 0.5);
    offsetUv -= floor(offsetUv);
    
    float height2 = perlin(offsetUv, rngState);
    
    vec2 weightXY = abs(uv - vec2(0.5)) * 2.0;
    float weight2 = weightXY.x * weightXY.y  + 0.001;
    float weight1 = (1.0-weightXY.x) * (1.0-weightXY.y) + 0.001;
    
    float weight = weight2 / (weight1 + weight2);
    weight = smoothstep(0., 1., weight);
    
    float v = mix(height1, height2, weight);*/
        

    
    // Increase contrast
    // (though this will flatten the tops and bottoms).
    v = v * 2.0 - 1.0;
    v = min(1.0, max(-1.0, v * 2.25));
    v = v * 0.5 + 0.5;
    
    // Make peaks pointer and valleys flatter.
    v *= v * v * v;
    
    v = 1. - v;
    
    
    fragColor = vec4(v, v, v, 1.0); 
    
    vec3 noise = getNoise(iFrame);
    fragColor = v > noise.z ? vec4(1) : vec4(0);
}
```


## Image

```glsl
// Controls:
// Mouse click - Clear and restart. (Useful for when you go fullscreen.)


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
float compress(float vec)
{
    float mag = abs(vec);

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
	vec2 uvC = fragCoord.xy / iChannelResolution[0].xy;
	vec2 uvR = (fragCoord.xy + vec2(1,0)) / iChannelResolution[0].xy;
	vec2 uvT = (fragCoord.xy + vec2(0,1)) / iChannelResolution[0].xy;
	vec2 uvL = (fragCoord.xy - vec2(1,0)) / iChannelResolution[0].xy;
	vec2 uvB = (fragCoord.xy - vec2(0,1)) / iChannelResolution[0].xy;
    vec2 colC = texture(iChannel0, uvC).zw;
    vec2 colR = texture(iChannel0, uvR).zw;
    vec2 colT = texture(iChannel0, uvT).zw;
    vec2 colL = texture(iChannel0, uvL).zw;
    vec2 colB = texture(iChannel0, uvB).zw;
    
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
    vec2 both = compress(gradient * 100.0 - curvature * 8000.0);
    
    
    // Unmodified gradients:
    // vec2 gradient = vec2(colC.r - colR.r, colC.r - colT.r);
    
    fragColor = vec4(both * 0.5 + 0.5, compress(colC.r), 1);
}

```
