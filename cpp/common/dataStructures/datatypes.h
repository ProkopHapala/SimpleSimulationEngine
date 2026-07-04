
#ifndef datatypes_h
#define datatypes_h

/// @file datatypes.h
/// @brief Bare POD structs for GPU buffers and C-style interchange — intentionally separate from Vec2T/Vec3T/Quat4T.
///
/// The Vec*T templates in math/ have methods, unions, and aliases — great for CPU code but
/// problematic for GPU kernel arguments and C-style file I/O where layout must be exactly
/// what you expect with no compiler-specific union padding. These structs (float4, int4, ...)
/// are guaranteed-layout POD types that map directly to OpenCL float4/int4 and to raw binary.
///
/// float8 and float16 exist because OpenCL wavefronts and SIMD groups naturally process
/// these widths — e.g. a single particle's full state (position, velocity, force, mass) fits
/// in a float8 or float16. double8 is the natural width for VertT<double> (8 doubles = 64 bytes).

template <typename T> struct vec2{
    union{
        struct{ T x,y; };
        T array[2];
    };
};

template <typename T> struct vec3{
    union{
        struct{ T x,y,z; };
        T array[3];
    };
};

template <typename T> struct vec4{
    union{
        struct{ T x,y,z,w;      };
        struct{ vec3<T> f; T e; };
        struct{ vec2<T> hi,lo;  };
        T array[4];
    };
};


//======= int

struct int2   { int x,y; };
struct int4   { int x,y,z,w; };
struct int8   { int x,y,z,w,hx,hy,hz,hw; };

struct double2 { double x,y; };
struct double4 { double x,y,z,w; };
struct double8 { double x,y,z,w,hx,hy,hz,hw; };

struct size_t4 { size_t x,y,z,w; };

//======= float

struct float2 {
	union{
		struct{ float x,y; };
		float array[2];
	};
};

struct float3 {
	union{
		struct{ float x,y,z; };
		float array[3];
	};
};

struct float4 {
	union{
		struct{ float x,y,z,w;     };
		struct{ float3 f; float e; };
		float array[4];
	};
};

struct float8 {
	union{
		struct{ float x,y,z,w,hx,hy,hz,hw; };
		struct{ float4 lo,hi; };
		float array[8];
	};
};

struct float16 {
	union{
		struct{ float x,y,z,w,hx,hy,hz,hw, x2,y2,z2,w2,hx2,hy2,hz2,hw2; };
		struct{ float8 lo,hi; };
		float array[16];
	};
};

#endif
