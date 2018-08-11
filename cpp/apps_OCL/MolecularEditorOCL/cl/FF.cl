
#define R2SAFE          1e-4f
#define RSAFE           1.0e-4f
#define COULOMB_CONST   14.399644f  // [eV/e]

float4 getCoulomb( float4 atom, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0f/( dot(dp,dp) +  R2SAFE );
     float   ir  = sqrt(ir2);
     float   E   = atom.w*sqrt(ir2);
     return (float4)(dp*(E*ir2), E );
}

float4 getLJ( float3 apos, float2 cLJ, float3 pos ){
     float3  dp  =  pos - apos;
     float   ir2 = 1.0f/( dot(dp,dp) + R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   E   =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  F   = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     return (float4)(F, E);
}

float4 getMorse( float3 dp, float3 REA ){
    //float3  dp  =  pos - apos;
    float   r     = sqrt( dot(dp,dp) + R2SAFE );
    float   expar = exp( REA.z*(r-REA.x) );
    float   E     = REA.y*expar*( expar - 2 );
    float   fr    = REA.y*expar*( expar - 1 )*2*REA.z;
    return (float4)(dp*(fr/r), E);
}

float8 getLJC( float4 atom, float2 cLJ, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0/( dot(dp,dp) +  R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   ELJ =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  FLJ = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     float   ir  = sqrt(ir2);
     float   Eel = atom.w*sqrt(ir2);
     return (float8)(FLJ, ELJ, dp*(Eel*ir2), Eel );
}

__kernel void evalCoulomb(
    int nAtoms, 
    __global float4* atoms, 
    __global float4*  poss,
    __global float4*    FE
){
    __local float4 LATOMS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        if(i>=nAtoms) break;
        //if(iL==0) printf("%i (%f,%f,%f)  %f \n", i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].w );
        
        LATOMS[iL] = atoms[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            fe += getCoulomb( LATOMS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
    //FE[iG] = poss[iG];
}

__kernel void evalLJ(
    const int nAtoms, 
    __global float4* atoms,
    __global float2*  cLJs,
    __global float4*  poss,
    __global float4*    FE
){
    __local float4 LATOMS[32];
    __local float2 LCLJS  [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        if(i>=nAtoms) break;
        //if(iL==0) printf("%i (%f,%f,%f)  %f \n", i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].w );
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            fe += getLJ( LATOMS[j].xyz, LCLJS [j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
    //FE[iG] = poss[iG];
}


__kernel void evalLJC(
    int nAtoms, 
    __global float4*   atoms,
    __global float2*    cLJs,
    __global float4*    poss,
    __global float8*    FE
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    fe.hi  = fe.hi*COULOMB_CONST;
    FE[iG] = fe;
}

__kernel void evalMorse(
    const int nAtoms, 
    __global float4*   atoms,
    __global float4*   REAs,
    __global float4*   poss,
    __global float4*   FE
){
    __local float4 lATOMs[32];
    __local float4 lREAs [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break; // wrong !!!!
        lATOMs[iL] = atoms[i];
        lREAs [iL] = REAs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getMorse( pos - lATOMs[j].xyz, lREAs[j].xyz );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    FE[iG] = fe;
}



__kernel void evalPLE(
    const int   nAtoms,
    const uint4 nGrid,
    const uint4 npbc,
    const float4 pos0,
    const float4 dA,
    const float4 dB,
    const float4 dC,
    const float alpha,
    __global float8*   atoms,
    __global float4*   FFPauli,
    __global float4*   FFLondon,
    __global float4*   FFelec 
){
    __local float8 lATOMs[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);

    // determine grid point position
    if(iG > (nGrid.x*nGrid.y*nGrid.z) ) return;
    uint ia =  iG %  nGrid.x;
    uint ib = (iG /  nGrid.x) % nGrid.y;
    uint ic =  iG / (nGrid.x*nGrid.y);
    float3 pos = pos + dA.xyz*ia + dB.xyz*ib + dC.xyz*ic;

    float3 A = dA.xyz*nGrid.x;
    float3 B = dB.xyz*nGrid.y;
    float3 C = dC.xyz*nGrid.z;

    float4 fp  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    float4 fl  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break; // wrong !!!!
        lATOMs[iL] = atoms[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){
                float3 dpos0 = pos - lATOMs[j].xyz;
                float3 REQi  = lATOMs[j].s456;
                // fe += getMorse( pos - lATOMs[j].xyz, lREAs[j].xyz );
                for (int jc=-npbc.z; jc<=npbc.z; jc++){
                    for (int jb=-npbc.z; jb<=npbc.y; jb++){
                        float3 dp = dpos0 + C*jc + B*jb + A*-npbc.z;
                        for (int ja=-npbc.z; ja<=npbc.x; ja++){

                            float   r    = sqrt( dot(dp,dp) );
                            float  ir    = 1.0f/(r+RSAFE);
                            float expar  = exp( alpha*(r-REQi.x) );
                            float fexp   = alpha*expar*REQi.y*ir;
                            fp.xyz       += dp * ( -fexp*expar*2.0f );                // repulsive part of Morse
                            fl.xyz       += dp * ( -fexp           );                // attractive part of Morse
                            fe.xyz       += dp * ( COULOMB_CONST*REQi.z*ir*ir*ir );  // Coulomb
                            // ToDo ... Energy as well ?
                            dp += A;
                        }
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FFPauli [iG] = fp;
    FFLondon[iG] = fl;
    FFelec  [iG] = fe;
}




