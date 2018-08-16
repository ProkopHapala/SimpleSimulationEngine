
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

    //FFPauli [iG] = (float4)(iG,iG,iG,iG);
    //FFLondon[iG] = (float4)(iG,iG,iG,iG);
    //FFelec  [iG] = (float4)(iG,iG,iG,iG);
    //if( (iG%100) == 0 ) printf("ig %i \n", iG);

    if( iG==0 ) printf("npbc (%i,%i,%i) \n", npbc.x, npbc.y, npbc.z );
    if( iG==0 ) printf("pos0 (%f,%f,%f) \n", pos0.x, pos0.y, pos0.z );
    if( iG==0 ) printf("dA (%f,%f,%f)   \n",   dA.x,   dA.y,   dA.z );
    if( iG==0 ) printf("dB (%f,%f,%f)   \n",   dB.x,   dB.y,   dB.z );
    if( iG==0 ) printf("dC (%f,%f,%f)   \n",   dC.x,   dC.y,   dC.z );


    const int n0x = -npbc.x;
    const int n1x =  npbc.x+1;
    const int n0y = -npbc.y;
    const int n1y =  npbc.y+1;
    const int n0z = -npbc.z;
    const int n1z =  npbc.z+1;
    
    /*
    if( iG==0 ){
        const int n0z = -npbc.z;
        const int n1z =  npbc.z;
        printf("jc %i .. %i \n", n0z, n1z );
        //for (int jc=0; jc<=npbc.z; jc++) printf("jc %i \n", jc );
        //for (int jc=-1; jc<=5; jc++) printf("jc %i \n", jc );
        //for (int jc=n0z; jc<=npbc.z; jc++) printf("jc %i \n", jc );
        for (int jc=n0z; jc<=n1z; jc++) printf("jc %i \n", jc );
        //for (int jc=-npbc.z; jc<=npbc.z; jc++) printf("jc %i \n", jc );
        //for (int jc=-npbc.z; jc<(npbc.z+1); jc++) printf("jc %i \n", jc );
        //for (int jc=-1; jc<(npbc.z*2+1); jc++) printf("jc %i \n", jc );
        return;
    }
    */

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break; // wrong !!!!
        lATOMs[iL] = atoms[i];
        //if( (iG%1000) == 0 ) printf("nAtoms %i i %i xyz (%f,%f,%f) req (%f,%f,%f)  \n", nAtoms, i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].s4, atoms[i].s5, atoms[i].s6 );
        if( iG == 0 ) printf("nAtoms %i i %i xyz (%f,%f,%f) req (%f,%f,%f)  \n", nAtoms, i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].s4, atoms[i].s5, atoms[i].s6 );
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( iG == 0 ) printf("i0 %i j %i i0+j %i \n", i0, j, i0+j );
            if( (j+i0)<nAtoms ){
                float3 dpos0 = pos - lATOMs[j].xyz;
                float3 REQi  = lATOMs[j].s456;

                if( iG == 0 ) printf("i dpos0 (%f,%f,%f) REQi (%f,%f,%f) \n", dpos0.x, dpos0.y, dpos0.z, REQi.x, REQi.y, REQi.z );
                if( iG == 0 ) printf("i %i (%i,%i,%i) .. (%i,%i,%i) \n", i0+j, n0x,n0y,n0z,   n1x,n1y,n1z );
                // fe += getMorse( pos - lATOMs[j].xyz, lREAs[j].xyz );
                //for (int jc=-npbc.z; jc<=npbc.z; jc++){ // for some reason does not work
                for (int jc=n0z; jc<n1z; jc++){
                    //if( iG == 0 ) printf("jc %i \n", jc );
                    //if( iG == 0 ) printf("jc %i \n", jc );
                    //for (int jb=-npbc.y; jb<=npbc.y; jb++){
                    for (int jb=n0y; jb<n1y; jb++){
                        float3 dp = dpos0 + C*jc + B*jb + A*n0x;
                        //if( iG == 0 ) printf("jbc (%i,%i) \n", jb,jc );
                        //for (int ja=-npbc.x; ja<=npbc.x; ja++){
                        for (int ja=n0x; ja<n1x; ja++){

                            //if( iG == 0 ) printf("ipbc (%i,%i,%i) \n", ja,jb,jc );
                            //if( iG == 0 ) printf("ipbc (%i,%i,%i) (%f,%f,%f) \n", ja,jb,jc, dp.x, dp.y, dp.z );

                            float   r    = sqrt( dot(dp,dp) );
                            float  ir    = 1.0f/(r+RSAFE);
                            float expar  = exp( alpha*(r-REQi.x) );
                            float fexp   = alpha*expar*REQi.y*ir;

                            if( iG==0 ) printf("ia %i jpbc (%i,%i,%i) r %f expar %e fexp %e \n", i0+j, ja,jb,jc,  r, expar, fexp );

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

    //if( (iG%100) == 0 ) printf("fp[%i] (%f,%f,%f,%f) \n", fp.x, fp.y, fp.z, fp.w);

    FFPauli [iG] = fp;
    FFLondon[iG] = fl;
    FFelec  [iG] = fe;

}




