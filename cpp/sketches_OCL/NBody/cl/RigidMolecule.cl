
#define R2SAFE 1.0e-2
#define F2MAX  10.0

#define R_MAX  1.8
#define R2MAX  3.24

float2 pair_force( float2 p1, float2 p2 ){
    float2 d   = p2 - p1;
    float  r2  = dot(d,d);
    if( r2 > R2MAX ) return (float2) (0.0f, 0.0f);
    float ir2 = 1/( r2 + R2SAFE );
    float fr  = (0.7-ir2)*(R2MAX-r2);
    return d * fr;
}

__kernel void NBody_force(
    unsigned int num,
    __global float2* pos, 
    __global float2* force
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local  float2 pos_shared[64];
    float2 p = pos[iG];
    float2 f = (float2) (0.0f, 0.0f);
    for (int i0=0; i0<num; i0+= nL ){ 
        pos_shared[iL] = pos[i0 + iL]; 
        barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
        for (int j=0; j<nL; j++){
            float2 pj  = pos_shared[j];
            f         += pair_force( p, pj );
        }
        barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
    }
    force[iG] = f;
}

float4 Quat_forceFromPoint( float4 q, float3 p, float3 fp ){
	
	float4 fq;
	float3 b;
	b.x =    p.y*q.y +               p.z*q.z;
	b.y =    p.x*q.y - 2*p.y*q.x +   p.z*q.w;
	b.z =    p.x*q.z -   p.y*q.w - 2*p.z*q.x;
	fq.x = dot(fp,b);

	b.x = -2*p.x*q.y +   p.y*q.x -   p.z*q.w;
	b.y =    p.x*q.x +               p.z*q.z;
	b.z =    p.x*q.w +   p.y*q.z - 2*p.z*q.y;
	fq.y = dot(fp,b);

	b.x = -2*p.x*q.z +   p.y*q.w +   p.z*q.x;
	b.y =   -p.x*q.w - 2*p.y*q.z +   p.z*q.y;
	b.z =    p.x*q.x +   p.y*q.y;
	fq.z = dot(fp,b);

	b.x =    p.y*q.z - p.z*q.y;
	b.y =   -p.x*q.z + p.z*q.x;
	b.z =    p.x*q.y - p.y*q.x;
	fq.w = dot(fp,b);
	return fq;
}

float4 Quat_qmul( float4 a, float4 b) {  // 16 mul, 12 add
    return (float4)(
      b.x * a.w + b.y * a.z - b.z * a.y + b.w * a.x ,
     -b.x * a.z + b.y * a.w + b.z * a.x + b.w * a.y ,
      b.x * a.y - b.y * a.x + b.z * a.w + b.w * a.z ,
     -b.x * a.x - b.y * a.y - b.z * a.z + b.w * a.w );
};

float3 Quat_transformVec( float4 q, float3 p){  
    return Quat_qmul( q, Quat_qmul( (float4)(q,0), q) );   // 32 mul, 24 add ... using matrix it would be faster (9 mul) (?)
}

float8 transformMol( int na, int nb, __local float4* atomsA, __local float4* atomsB ){
    float8 f;
    for(int i=0; i<na; i++){
        float3 pi = atomsA[i];
        float3 fi = 0.0f;
        for(int j=0; j<na; j++){
             fi += pair_force( pi, atomsB[j] );
        }
        Quat_forceFromPoint( q, pi, fp );
    }
    return f;
}

float8 transformMol( float8 trans, int n,  __local float4* fromAtom, __local float4* toAtoms ){
    for(int i=0; i<n; i++){ toAtoms[i] = Quat_transformVec( trans.hi, fromAtom[i] ) + trans.xyz; }
}

__kernel void RigidMolForce(
    int nAtoms,
    __global float8* pos,    // pos,quat 
    __global float8* force,  
    __const float4 AtomsInMol 
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local           Lpos[16];
    __local  float4 Latoms[16];
    __local  float4 atomsA[16];
    private  float4 atomsB[16]; // do we have enought memory for this ?
    Latoms[iL] = AtomsInMol[iL];
    barrier(CLK_LOCAL_MEM_FENCE);
    float8 p = pos[iG];
    float8 f = 0.0f;
    transformMol( p, nAtoms, Latoms, atomsB);
    for (int i0=0; i0<get_local_size(0); i0+= nL ){ 
        Lpos[iL] = pos[i0 + iL];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            transformMolMol( Lpos[j], nAtoms, Latoms, atomsB);
            f    += forceFromMol( nAtoms, nAtoms, atomsA, atomsB );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    force[iG] = f;
}



