
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

vec4 cubeFetch(samplerCube cube, ivec2 coord, int mip)
{
    vec2 uv = (vec2(coord) + vec2(0.5)) / vec2(textureSize(cube, mip));
    uv -= floor(uv);
    uv = uv * 2.0 - 1.0;
    
    vec3 ray = vec3(uv, 1);
    
    return textureLod(cube, ray, float(mip));
}

vec4 cubeLod(samplerCube cube, vec2 uv, float mip)
{
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



///////////////////////////////////////////////
// DATA
///////////////////////////////////////////////

#define DATA iChannel3

const int DataCompleteFrame = 4;

// Sampling Order -------------------------------
// Frame 0: Imprecise version calculated.
// Frane 2: Precise version calculated.
//     Depends on:
//         - Sampling Order, Frame 0
//         - Sampling Distance, Frame 1
// Frame 3: Ready.

const int SamplingOrderRow = 0;
const int SamplingOrderCount = 101;
const int SamplingOrderFrameImprecise = 0;
const int SamplingOrderFramePrecise = 2;


ivec3 GetSamplingOrderCoord(int n, sampler2D data)
{
    ivec2 coord = ivec2(n, SamplingOrderRow);
    
    return ivec3(texelFetch(data, coord, 0));
}

vec4 DebugSamplingOrderCoord(int n, sampler2D data)
{
    ivec3 d = GetSamplingOrderCoord(n, data);
    
    
    vec4 color = vec4
    (
        vec3(d.xyz) / vec3(vec2(10), 5), 
        1
    );
    
    color.xyz /= (1.0 + abs(color.xyz));
    color.xy += 0.5;
    return color;
}

// Cached Sampling Distances -------------------------------
// Frame 1: Calculate minimum distance between Pixel[0] and Pixel[N].
//     Depends on:
//         - Sampling Order, Frame 0
// Frame 3: Recalculate distances with new precise ordering.
//     Depends on:
//         - Sampling Order, Frame 2
const int SamplingDistanceRow = 1;
const int SamplingDistanceFrameImprecise = 1;
const int SamplingDistanceFramePrecise = 3;
float GetSamplingDistance(int n, sampler2D data)
{
    ivec2 coord = ivec2(n, SamplingDistanceRow);
    
    return texelFetch(data, coord, 0).x;
}

vec4 DebugSamplingDistance(int n, sampler2D data)
{
    float d = GetSamplingDistance(n, data);
    
    d *= 0.01;
    d = d / (1.0 + d);
    
    vec4 color = vec4(1);
    color.rgb = hsv2rgb(vec3(d, 1.0, 1.0));
    return color;
}



void debugAllData( sampler2D dataChannel, out vec4 fragColor, in vec2 fragCoord )
{      
    ivec2 iCoord = floorToInt(fragCoord);
    
    switch(iCoord.y)
    {
        case SamplingOrderRow: fragColor = DebugSamplingOrderCoord(iCoord.x, dataChannel); return;
        case SamplingDistanceRow: fragColor = DebugSamplingDistance(iCoord.x, dataChannel); return;
    }
    
    fragColor = texelFetch(dataChannel, iCoord, 0);
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
// This is a really stupid work-around to enforcing power-of-two buffers... but whatever.

void mainCubemap( out vec4 fragColor, in vec2 fragCoord, in vec3 rayOri, in vec3 rayDir )
{
    vec2 uv = fragCoord / iChannelResolution[0].xy; 
    
    float thresh = sin(iTime) * 0.25 + 0.5;
    fragColor = texture(iChannel0, uv).r > thresh ? vec4(1) : vec4(0);
}
```

## Image

```glsl
// HOLY ---- I DID IT
// Broad crawl outward, then a recursive mipmap quadtree search.
// FAST runtime, FAST compilation (comparitavely).
// Should compile in about ~7 seconds, ~70 FPS @ 1000x563

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
    fragColor = vec4(0, 0, 0, 1);
    
    ivec2 maxMipSize = textureSize(iChannel0, maxMip); 
    
    // Set up mouse coordinates.
    vec2 mouseUV = fragCoord / iResolution.xx * 2.0;
    vec2 mouseTestCoord = mouseUV * vec2(maxMipSize);
    ivec2 mouseTestIntCoord = floorToInt(mouseTestCoord);
    bool mouseIsOccupied = isOccupied(mouseUV); 
    
    // Set up fragment coordinates.
    vec2 fragUV = fragCoord / iResolution.xx * 2.0;
    vec2 fragTestCoord = fragUV * vec2(maxMipSize);
    ivec2 fragTestIntCoord = floorToInt(fragTestCoord);
    
    // Color the fragment.
    {
        float density = 0.;
        int count = 0;
        for(int i = maxMip-RECURSION_DEPTH; i <= maxMip; i++)
        {
            ivec2 mipSize = textureSize(iChannel0, i);
            bool fragIsOccupied = containsEdge(false, ivec2(fragUV * vec2(mipSize)), i);   
            bool fragIsOccupied2 = containsEdge(true, ivec2(fragUV * vec2(mipSize)), i);   

            density += fragIsOccupied ? 1. : 0.;
            density += !fragIsOccupied2 ? 1. : 0.;
            count+=2;
            
            // Remove this to show all mip levels explored.
            break;
        }
            
        // Compress the range to 0-1
        density /= 10. + density;
        
        density /= float(count) / (float(10 + count)); 
        //density = smoothstep(0., 1., density);
        //density *= density;
        
        fragColor.rgb = hsv2rgb(vec3(density*density*density*density*2.0 + 0.5, 1.-density * 0.5, density* 0.5));
    }
    
    vec3 diff_distSqr = vec3(4.0, 0.0, 16.0);

    int iCoordOffset = GetSampleIndexOffset((vec2(mouseTestIntCoord) + vec2(0.5) - mouseTestCoord));
    for(int iCoord = -1; iCoord < 64; iCoord++)
    {
        ivec2 samplingCoord = mouseTestIntCoord;// + iCoord == -1 ? SampleCoords[iCoord + iCoordOffset];
        if(iCoord >= 0)
        {
            samplingCoord += SampleCoords[iCoord + iCoordOffset];
        }
        
        vec3 sampleDiffDistSqr = get_diff_minDistSqr(mouseTestCoord, samplingCoord);
        
        // Need to add a bias for some reason.
        if(sampleDiffDistSqr.z >= diff_distSqr.z + 1.5)
            break;
            
        if(sampleDiffDistSqr.z >= diff_distSqr.z)
            continue; 
            
        if(samplingCoord == fragTestIntCoord)
            fragColor += vec4(0.1, 0, 0, 0);
            
        if(!containsEdge(mouseIsOccupied, samplingCoord, maxMip))
            continue;
            
        /*if(sampleDistSqr < distSqrToSignChange)
        {
            distSqrToSignChange = sampleDistSqr;
        }*/
            
        SubSampleArgs a;
        a.DIFFDISTSQR = diff_distSqr;
        a.SAMPLECOORD = samplingCoord;
        a.MIP = maxMip;
        a.TESTCOORD = mouseTestCoord;
        IncreaseDepth(a);
        diff_distSqr = SubSample
        (
            a,
            mouseIsOccupied
        ); 
    }
         
    vec3 diff_dist = vec3(diff_distSqr.xy, sqrt(diff_distSqr.z) * (mouseIsOccupied ? -1.0 : 1.0));     
         
    drawSDF
    (
        diff_dist,
        fragColor
    );
    
    // Compress the range to 0-1
    //fragColor /= vec4(0.25) + fragColor;
    //fragColor = smoothstep(vec4(0), vec4(1), fragColor);
}

```
