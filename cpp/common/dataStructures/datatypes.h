
#ifndef datatypes_h
#define datatypes_h

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
