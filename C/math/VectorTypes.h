#ifndef  VectorTypes_h
#define  VectorTypes_h

#define typedef_Vec2T(name,T)       typedef union name{ struct{T x,y;     }; struct{T a,b;     }; struct{T i,j;     }; T arr[2]; } name;
#define typedef_Vec3T(name,T)       typedef union name{ struct{T x,y,z;   }; struct{T a,b,c;   }; struct{T i,j,k;   }; T arr[3]; } name;
#define typedef_Vec4T(name,T)       typedef union name{ struct{T x,y,z,w; }; struct{T a,b,c,d; }; struct{T i,j,k,l; }; T arr[4]; } name;

#define typedef_Mat2T(name,T,V)     typedef union name{ struct{V a,b;     }; T arr[4]; struct{T ax,ay,bx,by; }; } name;
#define typedef_Mat3T(name,T,V)     typedef union name{ struct{V a,b,c;   }; T arr[9]; struct{T ax,ay,az, bx,by,bz, cx,cy,cz; }; struct{T xx,xy,xz, yx,yy,yz, zx,zy,zz; }; } name;
#define typedef_Mat4T(name,T,V)     typedef union name{ struct{V a,b,c,d; }; T arr[16]; } name;

//#define T int;
//typedef union { struct{T x,y,z; }; struct{T a,b,c; }; struct{T i,j,k; }; T arr[3]; } name;

typedef_Vec2T(int2   ,int   );
typedef_Vec2T(float2 ,float );
typedef_Vec2T(double2,double);

typedef_Vec3T(int3   ,int   );
typedef_Vec3T(float3 ,float );
typedef_Vec3T(double3,double);

typedef_Vec4T(int4   ,int   );
typedef_Vec4T(float4 ,float );
typedef_Vec4T(double4,double);

//typedef_Mat2T(int3   ,int   );
typedef_Mat2T(float2x2 ,float ,float2);
typedef_Mat2T(double2x2,double,double2);

//typedef_Mat2T(int3   ,int   );
typedef_Mat3T(float3x3 ,float ,float3);
typedef_Mat3T(double3x3,double,double3);

//typedef_Mat2T(int3   ,int   );
typedef_Mat4T(float4x4 ,float ,float4);
typedef_Mat4T(double4x4,double,double4);

#endif