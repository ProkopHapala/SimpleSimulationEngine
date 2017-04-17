
#define R2SAFE   1.0e-8
#define F2MAX   10.0

#define R_MAX  1.8
#define R2MAX  3.24

float3  atomicForceSR( float3 dp ){
    float r2 = dot(dp,dp);
    if (r2 > R2MAX) return (float3)(0.0f,0.0f,0.0f);
    float ir2 = 1/( r2 + R2SAFE );
    float fr  = (0.7-ir2)*(R2MAX-r2);
    return dp * fr ;
}

float3 atomicForceLJ( float3 dp, float C6, float C12 ){
    float ir2 = 1/( dot(dp,dp) + R2SAFE );
    float ir6 = ir2*ir2*ir2;
    float fr  = (6*C6*ir6 - 12*C12*ir6*ir6)*ir2;
    return dp * fr ;
}

float3 atomicForceR24( float3 dp, float C2, float C4 ){
    float ir2 = 1/( dot(dp,dp) + R2SAFE );
    float fr  = C2*ir2 - C4*ir2*ir2;
    return  dp * fr ;
}

float3 atomicForceCoulomb( float3 dp, float kQQ ){
    float ir2 = 1/( dot(dp,dp) + R2SAFE );
    float ir  = sqrt(ir2); 
    float fr  = ir*ir2*kQQ; 
    return dp * fr;
}

float3 atomicForceSpring( float3 dp, float k ){
    return dp * k;
}


__kernel void NBody_force(
    unsigned int num,
    __global float4* pos, 
    __global float4* force
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local  float4 pos_shared[64];
    float3 p = pos[iG].xyz;
    float3 f = (float3) (0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<num; i0+= nL ){ 
        pos_shared[iL] = pos[i0 + iL]; 
        barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
        for (int j=0; j<nL; j++){
            float3 pj  = pos_shared[j].xyz;
            f         += atomicForceR24( pj-p, 1.0f, 3.0f );
        }
        barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
    }
    force[iG] = (float4)(f,0.0f);
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
	return fq*2;
}

float4 Quat_qmul( float4 a, float4 b) {  // 16 mul, 12 add
    return (float4)(
      b.x * a.w + b.y * a.z - b.z * a.y + b.w * a.x ,
     -b.x * a.z + b.y * a.w + b.z * a.x + b.w * a.y ,
      b.x * a.y - b.y * a.x + b.z * a.w + b.w * a.z ,
     -b.x * a.x - b.y * a.y - b.z * a.z + b.w * a.w );
};

//float3 Quat_transformVec( float4 q, float3 p){  
//    return Quat_qmul( q, Quat_qmul( (float4)(q,0), q) );   // 32 mul, 24 add ... using matrix it would be faster (9 mul) (?)
//}

float3 Quat_transformVec  ( float4 q, float3 v){
    return v + 2 * cross(cross(v, q.xyz) + q.w*v, q.xyz);
}

float3 Quat_untransformVec( float4 q, float3 v){
    float3 tmp = 2 * cross(q.xyz, v );
    return v + q.w * tmp + cross(q.xyz, tmp );
}

float8 molForce( int na, int nb, float4 qrot, __local float4* atoms0, __local float4* atomsA, __local float4* atomsB ){
    float8 f;
    for(int i=0; i<na; i++){
        float3 pi = atomsA[i].xyz;
        float3 fi = (float3) (0.0f, 0.0f, 0.0f);
        for(int j=0; j<nb; j++){
             //fi += atomicForceSR( atomsB[j].xyz - pi );
             fi += atomicForceR24( atomsB[j].xyz - pi, 1.0f, 3.0f );
        }
        //f += (float8)( fi, 0.0f, Quat_forceFromPoint( qrot, atoms0[i].xyz, fi ) );
        f.xyz += fi;
        f.hi  += Quat_forceFromPoint( qrot, atoms0[i].xyz, fi );
    }
    return f;
}

void transformMol( float8 trans, int n,  __local float4* fromAtom, __local float4* toAtoms ){
    for(int i=0; i<n; i++){ 
        float4 p = fromAtom[i];
        toAtoms[i] = (float4)( Quat_transformVec( trans.hi, p.xyz ) + trans.xyz, p.w ); 
    }
}

__kernel void RigidMolTransform(
    int nAtoms,
    __constant float4* AtomsInMol,
    __global   float8* pos,    
    __global   float4* atomsT
){
    __local  float4 atoms0[16];
    
    const int i  = get_global_id (0);
    const int iG = get_group_id  (0); // should be equal to number of molecules
    const int iL = get_local_id  (0);
    //const int nL = get_local_size(0); // should be equal nAtoms
    
    atoms0[iL] = AtomsInMol[iL];
    barrier(CLK_LOCAL_MEM_FENCE);
    float8 posi = pos   [iG];
    float4 atom = atoms0[iL];
    //printf( "%i (%g,%g,%g,%g) (%g,%g,%g)\n", i, posi.hi.x, posi.hi.y, posi.hi.z, posi.hi.w,   atom.x, atom.y, atom.z );
    atomsT[i] = (float4)( Quat_transformVec( posi.hi, atom.xyz ) + posi.xyz, atom.w );    
}

__kernel void RigidMolAtomForce(
    int nAtoms,
    __constant float4* AtomsInMol,
    __global   float8* pos,    
    __global   float4* fatomsT
){
    __local  float4 atoms0[16];
    __local  float4 atomsB[16];
    
    const int iL = get_local_id  (0);
    
    atoms0[iL] = AtomsInMol[iL];
    barrier(CLK_LOCAL_MEM_FENCE);
    float8 posi   = pos   [get_group_id  (0)];
    float3 atom0  = atoms0[iL].xyz;
    float3 atomi  = Quat_transformVec( posi.hi, atom0 ) + posi.xyz;
    float3 fatomi = 0.0f;
    
    for (int jmol=0; jmol<get_num_groups(0); jmol++ ){
        float8 posj   = pos[jmol];
        float4 atomj0 = atoms0[iL];
        atomsB[iL]    = (float4)( Quat_transformVec( posj.hi, atomj0.xyz ) + posj.xyz, atomj0.w );
        barrier(CLK_LOCAL_MEM_FENCE);
        if(jmol!=get_group_id(0)){
            for (int j=0; j<get_local_size(0); j++){
                //fatomi += atomicForceSR( atomsB[j].xyz - atomi );
                fatomi += atomicForceR24( atomsB[j].xyz - atomi, 1.0f, 3.0f );
                //if((iL==0)&&(get_group_id(0)==0) ){
                //    //printf( "(%i,%i) (%f,%f,%f) (%f,%f,%f) (%f,%f,%f) \n", jmol,j, atomi.x, atomi.y, atomi.z,  atomsB[j].x, atomsB[j].y, atomsB[j].z,  fatomi.x, fatomi.y, fatomi.z );
                //    printf( "(%i,%i) (%f,%f,%f) (%f,%f,%f) \n", jmol,j, atomsB[j].x, atomsB[j].y, atomsB[j].z,  fatomi.x, fatomi.y, fatomi.z );
                //}
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    fatomsT[ get_global_id (0) ] = (float4)( fatomi, 0.0f );
}


__kernel void RigidMolForce(
    int nAtoms,
    __constant float4* AtomsInMol,
    __global   float8* pos,    // pos,quat 
    __global   float8* force 
){
    __local  float4 atoms0[16];
    __local  float4 atomsB[16];
    
    const int iL = get_local_id  (0);
    
    atoms0[iL] = AtomsInMol[iL];
    barrier(CLK_LOCAL_MEM_FENCE);
    float8 posi   = pos   [get_group_id  (0)];
    float3 atom0  = atoms0[iL].xyz;
    float3 atomi  = Quat_transformVec( posi.hi, atom0 ) + posi.xyz;
    float3 fatomi = 0.0f;
  
    for (int jmol=0; jmol<get_num_groups(0); jmol++ ){
        float8 posj   = pos[jmol];
        float4 atomj0 = atoms0[iL];
        atomsB[iL]    = (float4)( Quat_transformVec( posj.hi, atomj0.xyz ) + posj.xyz, atomj0.w );
        barrier(CLK_LOCAL_MEM_FENCE);
        if(jmol!=get_group_id(0)){
            for (int j=0; j<get_local_size(0); j++){
                //fatomi += atomicForceSR( atomsB[j].xyz - atomi );
                fatomi += atomicForceR24( atomsB[j].xyz - atomi, 1.0f, 3.0f );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    // sum force on molecule
    atoms0[iL] = (float4)(fatomi,0.0f);  // re-use local buffers 
    atomsB[iL] = Quat_forceFromPoint( posi.hi, atom0.xyz, fatomi ); 
    //atoms0[iL] = 1.0f;  
    barrier(CLK_LOCAL_MEM_FENCE);
    if(iL==0){
        float8 fmol  = 0.0f;
        for (int j=0; j<get_local_size(0); j++){
            fmol.lo += atoms0[j];
            fmol.hi += atomsB[j];
            //fmol.hi  = j;
            //fmol.hi  = get_local_size(0);
        }
        force[get_group_id(0)] = fmol;
        //force[get_group_id(0)] = 1515.0f;
    };
}


// this kernell is probably not feasibl due to limitations of local memory
__kernel void RigidMolForce_bak(
    int nAtoms,
    __constant float4* AtomsInMol,
    __global   float8* pos,    // pos,quat 
    __global   float8* force 
    
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local  float8   Lpos[16];
    __local  float4 atoms0[16];
    __local  float4 atomsA[16];
    __local  float4 atomsB[16];
    
    atoms0[iL] = AtomsInMol[iL];
    
    barrier(CLK_LOCAL_MEM_FENCE);
    float8 posi = pos[iG];
    float8 f    = 0.0f;
    transformMol( posi, nAtoms, atoms0, atomsA);      // each work_item should transform one atom, not whole molecule since we do not have enough local memory 
    
    for (int i0=0; i0<get_global_size(0); i0+= nL ){  // may be problem if number of molecules not divisible by "nL"
        Lpos[iL] = pos[i0 + iL];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            transformMol( Lpos[j], nAtoms, atoms0, atomsB);
            f    += molForce( nAtoms, nAtoms, posi.hi, atoms0, atomsA, atomsB );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    force[iG] = f;
    
}



