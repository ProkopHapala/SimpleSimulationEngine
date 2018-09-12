
// https://www.khronos.org/registry/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html

#define R2ELEC 1.0
#define R2SAFE          1e-4f
#define RSAFE           1.0e-4f
#define COULOMB_CONST   14.399644f  // [eV/e]


#define PRECISE_INTERPOLATION

//__constant sampler_t sampler_1 = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_1 = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_2 = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

/*

If the sampler is specified as using unnormalized coordinates (floating-point or integer
coordinates), filter mode set to CLK_FILTER_NEAREST and addressing mode set to one of the
following modes - CLK_ADDRESS_NONE, CLK_ADDRESS_CLAMP_TO_EDGE or
CLK_ADDRESS_CLAMP, the location of the image element in the image given by (i, j, k) in
section 8.2 will be computed without any loss of precision.

For all other sampler combinations of normalized or unnormalized coordinates, filter and
addressing modes, the relative error or precision of the addressing mode calculations and the
image filter operation are not defined by this revision of the OpenCL specification. To ensure a
minimum precision of image addressing and filter calculations across any OpenCL device, for
these sampler combinations, developers should unnormalize the image coordinate in the kernel
and implement the linear filter in the kernel with appropriate calls to read_image{f|i|ui} with a
sampler that uses unnormalized coordinates, filter mode set to CLK_FILTER_NEAREST,
addressing mode set to CLK_ADDRESS_NONE, CLK_ADDRESS_CLAMP_TO_EDGE or
CLK_ADDRESS_CLAMP and finally performing the interpolation of color values read from the
image to generate the filtered color value. 

*/


inline float4 read_imagef_trilin( __read_only image3d_t imgIn, float4 coord ){
    float4 icoord;
    float4 fc     =  fract( coord, &icoord );
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    return  
     (( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(1.0,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,1.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(1.0,1.0,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,1.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(1.0,0.0,1.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,1.0,1.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(1.0,1.0,1.0,0.0) ) * fc.x )*fc.y )*fc.z;
}; 

inline float3 rotMat ( float3 v, float3 a, float3 b, float3 c ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }
inline float3 rotMatT( float3 v,  float3 a, float3 b, float3 c  ){ return a*v.x + b*v.y + c*v.z; }


inline float3 rotQuat( float4 q, float3 v ){
    // http://www.geeks3d.com/20141201/how-to-rotate-a-vertex-by-a-quaternion-in-glsl/
    //return v + 2.0 * cross(q.xyz, cross(q.xyz, v) + q.w * v);
    //https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    float3 t = 2 * cross( q.xyz, v );
    return v + t * q.w + cross( q.xyz , t );
    // https://gamedev.stackexchange.com/questions/28395/rotating-vector3-by-a-quaternion
    //return 2.0 * dot(q.xyz, v) * q.xyz + ( q.w*q.w - dot(q.xyz,q.xyz) ) * v + 2.0 * q.w * cross(q.xyz, v);
}

inline float4 drotQuat_exact( float4 q, float4 tq ){
    float r = length( tq );
    float a = r * 0.5;
    //float2 csa = (float2)(cos(a),sin(a)/r);
    tq.xyz *= sin(a)/(r+1.0e-8f);
    return q*cos(a) + (float4){
        tq.x*q.w + tq.y*q.z - tq.z*q.y,
       -tq.x*q.z + tq.y*q.w + tq.z*q.x,
        tq.x*q.y - tq.y*q.x + tq.z*q.w,
       -tq.x*q.x - tq.y*q.y - tq.z*q.z,
       // dot(tq.xy,q.wz) - tq.z*q.y,
       // dot(tq.yz,q.wx) - tq.x*q.z,
       // dot(tq.xz,q.yw) - tq.y*q.z,
       //-dot(tq,q.xyz),
    };
    //return q.wzyz*tq.x + q.zwxy*tq.y + q.yxwz*tq.z + q.xyzw*tq.w;
}

float3 tipForce( float3 dpos, float4 stiffness, float4 dpos0 ){
    float r = sqrt( dot( dpos,dpos) );
    return  (dpos-dpos0.xyz) * stiffness.xyz              // harmonic 3D
         + dpos * ( stiffness.w * (r-dpos0.w)/r );  // radial
}

float4 interpFE( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef( imgIn, sampler_1, coord );
}

// this should be macro, to pass values by reference
void move_LeapFrog( float3 f, float3 p, float3 v, float2 RP ){
    v  =  f * RP.x + v*RP.y;
    p +=  v * RP.x;
}

#define FTDEC 0.5f
#define FTINC 1.1f
#define FDAMP 0.99f
#define N_RELAX_STEP_MAX  64
#define F2CONV  1e-8

// this should be macro, to pass values by reference
void move_FIRE( float3 f, float3 p, float3 v, float2 RP, float4 RP0 ){
    // RP0 = (t?,damp0,tmin,tmax)
    float ff   = dot(f,f);
    float vv   = dot(v,v);
    float vf   = dot(v,f);
    if( vf < 0 ){
        v      = 0.0f;
        RP.x   = max( RP.x*FTDEC, RP0.z );    // dt
        RP.y   = RP0.y;                      // damp
    }else{
        v      = v*(1-RP.y) + f*RP.y * sqrt(vv/ff);
        RP.x   = min( RP.x*FTINC, RP0.w );   // dt
        RP.y  *= FDAMP;                     // damp
    }
    // normal leap-frog times step
    v +=  f * RP.x;
    p +=  v * RP.x;
}

__kernel void getFEgrid(
    __read_only image3d_t  imgPauli,
    __read_only image3d_t  imgLondon,
    __read_only image3d_t  imgElec,
    __global  float4*  PLQs,    // const memory ?
    __global  float4*  poss,
    __global  float4*  FEs,
    float4 pos0,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC
){
    const float3 pos   = poss[ get_global_id (0) ].xyz + pos0.xyz;
    const float3 PLQ   = PLQs[ get_global_id (0) ].xyz;
    float4 coord = (float4)( dot(pos,dinvA.xyz),dot(pos,dinvB.xyz),dot(pos,dinvC.xyz), 0.0f );
    float4 fe;
    #ifdef PRECISE_INTERPOLATION
    coord += (float4)(-0.5,-0.5,-0.5,0.0);
    fe  = PLQ.x * read_imagef_trilin( imgPauli,  coord );
    fe += PLQ.y * read_imagef_trilin( imgLondon, coord );
    fe += PLQ.z * read_imagef_trilin( imgElec,   coord );
    #else
    fe  = PLQ.x * read_imagef( imgPauli,  sampler_1, coord );
    fe += PLQ.y * read_imagef( imgLondon, sampler_1, coord );
    fe += PLQ.z * read_imagef( imgElec,   sampler_1, coord );
    #endif
    FEs[ get_global_id (0) ] = fe;

    //if( get_global_id (0)==0 ) printf( "GPU PLQ %g %g %g \n", PLQ.x, PLQ.y, PLQ.z );
    //printf( "GPU %i pos(%g,%g,%g) PLQ(%g,%g,%g) fe(%g,%g,%g) \n", get_global_id(0), pos.x, pos.y, pos.z, PLQ.x, PLQ.y, PLQ.z, fe.x,fe.y,fe.z );
}

/*
Purpose,
    have N systems each composed of M molecules each composed of K atoms
    - assume M,K < local_size (~32)
    They are laying on grid
    we want to relax all the systems in paralel
    later we also want to find which configurations overlaps with previous
*/

__kernel void getForceRigidSystemSurfGrid(
    __read_only image3d_t  imgPauli,
    __read_only image3d_t  imgLondon,
    __read_only image3d_t  imgElec,
    // found - previously found molecular configurations
    __global  int2*    mol2atoms,    // mol2atoms[type.x:type.y]
    //__global  int2*  confs,        // pointer to poses ... since we have constant number of molecules, we dont need this
    __global  float8*  atomsInTypes, // atoms in molecule types
    __global  float4*  PLQinTypes,   // atoms in molecule types
    __global  float8*  poses,        // pos, qrot
    __global  float8*  fposes,       // force acting on pos, qrot
    __global  float8*  vposes,       // force acting on pos, qrot
    float4 pos0,
    float4 dinvA,
    float4 dinvB,
    float4 dinvC,
    int nSystems,
    int nMols, // nMols should be approx local size
    float alpha,
    float dt,
    float damp,
    int nStep
){

    __local float8 lATOMs[32]; // here we store atoms of molecule j
    __local float8 lPOSES[8];  // here we store atoms of molecule j

    //const int iG = get_global_id (0);
    const int nL    = get_local_size(0);
    // we asume there is nMol and nAtoms < local_size
    const int imol  = get_local_id(0);
    const int iL    = get_local_id(0);
    const int isys  = get_group_id(0);
    
/*
    if( get_global_id(0)==0 ){
        printf( "GPU nGlobal %i nLocal %i nSystems %i nMols %i \n", (int)get_global_size(0), nL, nSystems, nMols );
        //printf( "dinvA (%g,%g,%g,%g)\n", dinvA.x, dinvA.y, dinvA.z );
        //printf( "dinvB (%g,%g,%g,%g)\n", dinvB.x, dinvB.y, dinvB.z );
        //printf( "dinvC (%g,%g,%g,%g)\n", dinvC.x, dinvC.y, dinvC.z);
        printf( "GPU nStep %i dt %f damp %f \n", nStep, dt, damp );
    }
*/

    //if( isys > nSystems ) return;  // perhaps not needed if global_size set properly
    //if( imol > nMols ) return;

    const int molOffset = isys * nMols;
    const int oimol     = molOffset+imol;
    const int iatomi    = mol2atoms[oimol].x;
    const int natomi    = mol2atoms[oimol].y;
    float4 mposi        = poses[oimol].lo;
    float4 qroti        = poses[oimol].hi;
    float4 forceE       = (float4)(0.0,0.0,0.0,0.0);
    float4 torq         = (float4)(0.0,0.0,0.0,0.0); // WHAT IS DIFFERENCE BETWEEN dqrot/dforce and torq?

    float4 forceV       = vposes[oimol].lo;
    float4 torqV        = vposes[oimol].hi; // WHAT IS DIFFERENCE BETWEEN dqrot/dforce and torq?

/*
    if( get_local_id(0)==0 ){
    //if( get_global_id(0)==0 ){ 
        //printf( "GPU isystem %i \n", isys );
        //printf( "GPU isys,oimol(%i,%i) i0,n(%i,%i) p(%g,%g,%g) q(%g,%g,%g,%g) \n", isys, oimol, iatomi, natomi,  mposi.x,mposi.y,mposi.z,  qroti.x,qroti.y,qroti.z,qroti.w  );
        printf( "GPU natomi %i \n", natomi );
        for(int ia=0; ia<natomi; ia++){
            float8 ai = atomsInTypes[iatomi+ia];
            printf( "GPU isys %i ia %i %i xyz(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", isys, imol, ia, iatomi+ia, ai.x, ai.y, ai.z, ai.w,  ai.s4, ai.s5, ai.s6, ai.s7 );
        }
    }
*/

    for(int iStep=0; iStep<nStep; iStep++ ){

        forceE = (float4)(0.0,0.0,0.0,0.0);
        torq   = (float4)(0.0,0.0,0.0,0.0);

        // ==== Molecule - Molecule Interaction

        lPOSES[iL].lo = mposi;
        lPOSES[iL].hi = qroti;
        barrier(CLK_LOCAL_MEM_FENCE);

        for(int jmol=0; jmol<nMols; jmol++){

            // transform atoms of jmol to world coords and store them local memory
            const int ojmol  = molOffset+jmol;
            const int iatomj = mol2atoms[ojmol].x;
            const int natomj = mol2atoms[ojmol].y;
            if(iL<natomj){
                float8 atomj = atomsInTypes[iatomj+iL];
                //atomj.xyz    = rotQuat( poses[ojmol].hi, atomj.xyz ) + poses[ojmol].xyz; // TODO: we store poses in local memory
                atomj.xyz    = rotQuat( lPOSES[jmol].hi, atomj.xyz ) + lPOSES[jmol].xyz; 
                lATOMs[iL]   = atomj;

                //if( (isys==1) )  printf( "GPU lATOM imol %i jmol %i ja %i pi %5.5e %5.5e %5.5e p %5.5e %5.5e %5.5e q %5.5e %5.5e %5.5e %5.5e \n", imol, jmol, iL, atomj.x,atomj.y,atomj.z, 
                //        poses[ojmol].x,poses[ojmol].y,poses[ojmol].z,
                //        poses[ojmol].hi.x,poses[ojmol].hi.y,poses[ojmol].hi.z,poses[ojmol].hi.w        
                //);

            }

            barrier(CLK_LOCAL_MEM_FENCE);
            if ( (jmol==imol) || (imol>=nMols) ) continue; // prevent self-interaction and reading outside buffer
            for(int ia=0; ia<natomi; ia++){ // atoms of molecule i
                float4 adposi = atomsInTypes[iatomi+ia].lo;
                adposi.xyz    = rotQuat( qroti, adposi.xyz );
                float3 aposi  = adposi.xyz + mposi.xyz;
                float3 REQi   = atomsInTypes[iatomi+ia].s456;

                // molecule-molecule interaction
                float4 fe = (float4)(0.0,0.0,0.0,0.0);

                for(int ja=0; ja<natomj; ja++){            // atoms of molecule j
                    float3 dp   = lATOMs[ja].xyz - aposi;  // already transformed
                    float3 REQj = lATOMs[ja].s456;
                    // force
                    float r0    = REQj.x + REQi.x;
                    float eps   = REQj.y * REQi.y; 
                    float cElec = REQj.z * REQi.z * COULOMB_CONST;
                    float r     = sqrt( dot(dp,dp) + R2SAFE );
                    float expar = exp( alpha*(r-r0));
                    float fr    = eps*2*alpha*( expar*expar - expar ) - cElec/( r*r + R2ELEC );
                    fe.xyz     += dp *( fr/r );
                    fe.w       += eps*( expar*expar - 2*expar ) + cElec/( r ); // Energy - TODO there should be approx of arctan(x)

                    //if( isys==0 ) printf( "GPU (%i,%i) (%i,%i) r %g expar %g fr %g kqq %g a %g eps %g \n", imol, ia, jmol, ja, r, expar, fr, cElec, alpha, eps );

                    //if( (isys==1) && (iStep==1) && (jmol==0) ) 
                    //if( (isys==1)  ) 
                    //if( (iStep==1) ) 
                    //if( (jmol==0) ) 
                    //if( (isys==1) && (jmol==0) ) 
                    //if (isys==1) printf( "GPU fm imol %i jmol %i ia %i ja %i dp  %5.5e %5.5e %5.5e pi %5.5e %5.5e %5.5e f %5.5e %5.5e %5.5e \n", imol, jmol, ia, ja, dp.x,dp.y,dp.z, adposi.x, adposi.y, adposi.z, fe.x,fe.y,fe.z );
                }

                forceE += fe;
                torq   += cross(adposi, fe);

            }

            //if(isys == 1 ) printf( "%i imol %i p  %5.5e %5.5e %5.5e f %5.5e %5.5e %5.5e \n", iStep, imol, mposi.x,mposi.y,mposi.z, forceE.x, forceE.y, forceE.z );

            barrier(CLK_LOCAL_MEM_FENCE);
        }

        if ( imol<nMols ){
            // ==== Molecule - Grid interaction

            // TODO: we may store transformed atoms of imol to local mem ?
            
            for(int ia=0; ia<natomi; ia++){ // atoms of molecule i
                const float4 adposi = atomsInTypes[iatomi+ia].lo;
                adposi.xyz    = rotQuat( qroti, adposi.xyz );
                const float3 aposi  = adposi.xyz + mposi.xyz + pos0.xyz;

                float4 fe = (float4)(0.0,0.0,0.0,0.0);
                const float4 coord = (float4)( dot(aposi,dinvA.xyz),dot(aposi,dinvB.xyz),dot(aposi,dinvC.xyz), 0.0f );

                float3 REQi    = atomsInTypes[iatomi+ia].s456;
                float expar    = exp(-alpha*REQi.x);
                float cPauli   =    REQi.y*expar*expar;
                float cLondon  = -2*REQi.y*expar;
                fe += cPauli *read_imagef( imgPauli,  sampler_1, coord );
                fe += cLondon*read_imagef( imgLondon, sampler_1, coord );
                fe += REQi.z *read_imagef( imgElec,   sampler_1, coord );
                
                //const float3 PLQi   = PLQinTypes[iatomi+ia].xyz;
                //fe += PLQi.x*read_imagef( imgPauli,  sampler_1, coord );
                //fe += PLQi.y*read_imagef( imgLondon, sampler_1, coord );
                //fe += PLQi.z*read_imagef( imgElec,   sampler_1, coord );
                //fe = read_imagef( imgLondon, sampler_1, coord );

                //if( isys==0 && imol==0 ) printf( "GPU fgrid: imol %i ia %i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, aposi.x, aposi.y, aposi.z,  fe.x, fe.y, fe.z );
                //if( isys==0 && imol==0 ) printf( "GPU fgrid: imol %i ia %i p(%5.5e,%5.5e,%5.5e) g(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, aposi.x, aposi.y, aposi.z,  coord.x,coord.y,coord.z,  fe.x, fe.y, fe.z );
                //if( isys==0 && imol==0 ) printf( "GPU fgrid: imol %i ia %i p(%5.5e,%5.5e,%5.5e) PLQ(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, aposi.x, aposi.y, aposi.z,  cPauli,cLondon,REQi.z,  fe.x, fe.y, fe.z );

                //if( isys==0 ) printf( "GPU fgrid: imol %i ia %i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, aposi.x, aposi.y, aposi.z,  fe.x, fe.y, fe.z );
                //if( isys==0 ) printf( "GPU fgrid: imol %i ia %i pg(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, coord.x, coord.y, coord.z,  fe.x, fe.y, fe.z );
                //if( (isys==0) && (imol==0) && (ia==0) ) printf( "GPU fgrid: iG %i imol %i ia %i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n",  
                //                                                 get_global_id(0), imol, ia, aposi.x, aposi.y, aposi.z,  fe.x, fe.y, fe.z);

                //if( (isys==0) && (imol==0) ) printf( "GPU fgrid: iG %i imol %i ia %i REPLQ(%5.5e,%5.5e,%5.5e) PLQ(%5.5e,%5.5e,%5.5e) \n",   get_global_id(0), imol, ia, cPauli, cLondon, REQi.z,  PLQi.x, PLQi.y, PLQi.z  );
                //if( (isys==0) && (imol==0) ) printf( "GPU imol %i ia %i RE(%5.5e,%5.5e)->PL(%5.5e,%5.5e) PL(%5.5e,%5.5e) \n", imol, ia, REQi.x, REQi.y, cPauli, cLondon, PLQi.x, PLQi.y  );

                forceE += fe;
                torq   += cross(adposi, fe);

                //if( isys==0 ) printf( "GPU imol %i f(%g,%g,%g) tq(%g,%g,%g) \n", imol, forceE.x, forceE.y, forceE.z, torq.x, torq.y, torq.z );
                //printf( "GPU isystem %i imol %i f(%g,%g,%g) tq(%g,%g,%g) \n", isys, imol, forceE.x, forceE.y, forceE.z, torq.x, torq.y, torq.z );

            } // ia

            // Gradient-descent step  moveSystemGD
            //dt = 0.1;
            if(dt>-1e-8){
                //if( isys==0 ) printf( "GPU fgrid: imol %i p(%5.5e,%5.5e,%5.5e) q(%5.5e,%5.5e,%5.5e,%5.5e) \n", imol, poses[oimol].x, poses[oimol].y, poses[oimol].z,  poses[oimol].s4, poses[oimol].s5, poses[oimol].s6, poses[oimol].s7 );
                //if( isys==0 ) printf( "GPU move: imol %i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, poses[oimol].x, poses[oimol].y, poses[oimol].z,  fposes[oimol].x, fposes[oimol].y, fposes[oimol].z );
                //if( isys==0 && imol==0 ) printf( "GPU move: step %i imol %i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", iStep, imol, mposi.x, mposi.y, mposi.z,  forceE.x, forceE.y, forceE.z );

                //if(isys == 1 && imol==1 ) printf( "%i p(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) v(%5.5e,%5.5e,%5.5e) \n", iStep, mposi.x,mposi.y,mposi.z, forceE.x, forceE.y, forceE.z, forceV.x, forceV.y, forceV.z );
                //if(isys == 1 && imol==1 ) printf( "%i p %5.5e %5.5e %5.5e f %5.5e %5.5e %5.5e v %5.5e %5.5e %5.5e \n", iStep, mposi.x,mposi.y,mposi.z, forceE.x, forceE.y, forceE.z, forceV.x, forceV.y, forceV.z );

                //if(isys == 1 ) printf( "%i imol %i p  %5.5e %5.5e %5.5e f %5.5e %5.5e %5.5e v %5.5e %5.5e %5.5e \n", iStep, imol, mposi.x,mposi.y,mposi.z, forceE.x, forceE.y, forceE.z, forceV.x, forceV.y, forceV.z );

                //if(isys == 1 && imol==1 ) printf( "%i |f| %5.5e |tq| %5.5e p(%5.5e,%5.5e,%5.5e)\n", iStep, dot( forceE.xyz,forceE.xyz), dot(torq.xyz,torq.xyz), mposi.x,mposi.y,mposi.z );
                //if(isys == 1 && imol==1 ) printf( "%i |f| %5.5e |tq| %5.5e p(%5.5e,%5.5e,%5.5e)\n", iStep, dot( forceE.xyz,forceE.xyz), dot(torq.xyz,torq.xyz), mposi.x,mposi.y,mposi.z );

                //mposi.xyz += forceE.xyz * dt;
                //qroti      = drotQuat_exact( qroti, torq * dt );

                forceV.xyz = forceV.xyz*damp + forceE.xyz * dt;
                torqV.xyz  = torqV.xyz *damp + torq.xyz   * dt;
                mposi.xyz += forceV.xyz * dt;
                qroti      = drotQuat_exact( qroti, torqV * dt );

            }
        } // imol<nMols

    } // iStep;

    if ( imol<nMols ){
        fposes[oimol].lo  = forceE;
        fposes[oimol].hi  = torq;
        vposes[oimol].lo  = forceV;
        vposes[oimol].hi  = torqV;
         poses[oimol].lo  = mposi;
         poses[oimol].hi  = qroti;
    }

}









