
#ifndef datatypes_h
#define datatypes_h

//======= int

struct int2   { int x,y; };

//======= float

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

#endif
