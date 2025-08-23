/*

// ============= For automatic generation of interfaces

int nnode
int natom
int nvec = nnode+natom

//------- Dynamical

_RW float4*  apos      [ nSys* nvec ]
_RW float4*  aforce    [ nSys* nvec ]     
_RW float4*  avel      [ nSys* nvec ]    
_RW float4*  fneigh    [ nSys* nnode*2*4 ]

//------- parameters

_R  int4*    neighs    [ nSys* natom ]    
_R  int4*    bkNeighs  [ nSys* natom ]    
_R  int4*    neighCell [ nSys* natom ]    
_R  float4*  REQs     [ nSys* natom ]
_R  float4*  apars,    [ nSys* nnode ] 
_R  float4*  bLs,      [ nSys* nnode ] 
_R  float4*  bKs,      [ nSys* nnode ] 
_R  float4*  Ksp,      [ nSys* nnode ]
_R  float4*  Kpp,      [ nSys* nnode ]
_R  cl_Mat3* lvecs,    [ nSys ]
_R  cl_Mat3* ilvecs    [ nSys ]

_R float4*  apos_surf   [ natom_surf ]
_R float4*  aforce_surf [ natom_surf ]  

_R float4*  ps,         [ ndipol ]    
_R float4*  dipols,     [ ndipol ]

_RW image3d_t  FE_Paul[ ng.x, ng.y, ng.z ]
_RW image3d_t  FE_Lond[ ng.x, ng.y, ng.z ]
_RW image3d_t  FE_Coul[ ng.x, ng.y, ng.z ]

const float4   MDpars
const int4     nGrid
const cl_Mat3  dGrid
const float4   grid_p0 

*/


/*
    relax_multi.cl -  OpenCL kernel source code for multi-system relaxation
    
    The file contains various OpenCL kernels icluding:

    
    This file is part of FireCode project. 

*/


#define iGdbg 0
#define iSdbg 0

#pragma OPENCL EXTENSION cl_khr_fp64 : disable
//#pragma OPENCL FP_CONTRACT ON

// ======================================================================
// ======================================================================
//                           FUNCTIONS
// ======================================================================
// ======================================================================

//#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

typedef struct __attribute__ ((packed)){
    float4 a;
    float4 b;
    float4 c;
} cl_Mat3;

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define const_kB        8.617333262145e-5f   // [ eV/K ]



inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2){  a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y }; }     // divison of unitary complex numbers (i.e. rotation backwards)
//inline void     udiv_cmplx(               const VEC& b ){                            T x_ =    x*b.x +   y*b.y;         y =    y*b.x -   x*b.y;       x=x_;  }

inline float3 rotMat ( float3 v, float3 a, float3 b, float3 c ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }  // rotate vector v by matrix (a,b,c)
inline float3 rotMatT( float3 v, float3 a, float3 b, float3 c ){ return a*v.x + b*v.y + c*v.z; }                 // rotate vector v by matrix (a,b,c) transposed

//<<<file Forces.cl

// evaluate angular force and energy using cos(angle) formulation,    - faster, but not good for angles > 90 deg
inline float evalAngCos( const float4 hr1, const float4 hr2, float K, float c0, __private float3* f1, __private float3* f2 ){
    float  c = dot(hr1.xyz,hr2.xyz);
    float3 hf1,hf2;
    hf1 = hr2.xyz - hr1.xyz*c;
    hf2 = hr1.xyz - hr2.xyz*c;
    float c_   = c-c0;
    float E    = K*c_*c_;
    float fang = -K*c_*2;
    hf1 *= fang*hr1.w;
    hf2 *= fang*hr2.w;
    *f1=hf1;
    *f2=hf2;
    return E;
}

// evaluate angular force and energy using cos(angle/2) formulation - a bit slower, but not good for angles > 90 deg
inline float evalAngleCosHalf( const float4 hr1, const float4 hr2, const float2 cs0, float k, __private float3* f1, __private float3* f2 ){
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    // the main advantage is that it is quasi-harmonic beyond angles > 90 deg
    float3 h  = hr1.xyz + hr2.xyz;  // h = a+b
    float  c2 = dot(h,h)*0.25f;     // cos(a/2) = |ha+hb|  (after normalization)
    float  s2 = 1.f-c2 + 1e-7f;      // sin(a/2) = sqrt(1-cos(a/2)^2) ;  s^2 must be positive (otherwise we get NaNs)
    float2 cso = (float2){ sqrt(c2), sqrt(s2) }; // cso = cos(a/2) + i*sin(a/2)
    float2 cs = udiv_cmplx( cs0, cso );          // rotate back by equilibrium angle
    float  E         =  k*( 1 - cs.x );          // E = k*( 1 - cos(a/2) )  ; Do we need Energy? Just for debugging ?
    float  fr        = -k*(     cs.y );          // fr = k*( sin(a/2) )     ; force magnitude
    c2 *= -2.f;
    fr /=  4.f*cso.x*cso.y;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a) 
    float  fr1    = fr*hr1.w; // magnitude of force on atom a
    float  fr2    = fr*hr2.w; // magnitude of force on atom b
    *f1 =  h*fr1  + hr1.xyz*(fr1*c2);  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| ); force on atom a
    *f2 =  h*fr2  + hr2.xyz*(fr2*c2);  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| ); force on atom b
    return E;
}

// evaluate angular force and energy for pi-pi alignment interaction
inline float evalPiAling( const float3 h1, const float3 h2,  float K, __private float3* f1, __private float3* f2 ){  // interaction between two pi-bonds
    float  c = dot(h1,h2); // cos(a) (assumes that h1 and h2 are normalized)
    float3 hf1,hf2;        // working forces or direction vectors
    hf1 = h2 - h1*c;       // component of h2 perpendicular to h1
    hf2 = h1 - h2*c;       // component of h1 perpendicular to h2
    bool sign = c<0; if(sign) c=-c; // if angle is > 90 deg we need to flip the sign of force
    float E    = -K*c;     // energy is -K*cos(a)
    float fang =  K;       // force magnitude
    if(sign)fang=-fang;    // flip the sign of force if angle is > 90 deg
    hf1 *= fang;           // force on atom a
    hf2 *= fang;           // force on atom b
    *f1=hf1;
    *f2=hf2;
    return E;
}

// evaluate bond force and energy for harmonic bond stretching
inline float evalBond( float3 h, float dl, float k, __private float3* f ){
    float fr = dl*k;   // force magnitude
    *f = h * fr;       // force on atom a
    return fr*dl*0.5f;  // energy
}

// evaluate non-covalent interaction force and energy for Lennard-Jones (Q) and Coulomb interactions of charges (Q) and hydrogen bond correction (pseudo-charges H), damping R2damp is used to avoid singularity at r=0
inline float4 invR2( float3 dp ){
    const float ir2 = 1.f/(dot(dp,dp));
    const float E   = ir2;
    return  (float4){ dp*ir2*ir2, E };
}

inline float4 R2gauss( float3 dp ){
    const float r2 = dot(dp,dp);
    if(r2>1.0){ return (float4){0.f,0.f,0.f,0.f}; }
    float p = 1 - r2;
    return  (float4){ dp*p, p*p }; // dp*r2 = |dp|
}


// evaluate non-covalent interaction force and energy for Lennard-Jones (Q) and Coulomb interactions of charges (Q) and hydrogen bond correction (pseudo-charges H), damping R2damp is used to avoid singularity at r=0
inline float4 exp_r( float3 dp, float b ){
    const float r = length(dp);
    const float E = exp(-b*r);
    return  (float4){ dp*(E*b/r), E };
}

inline float4 exp_r_lin4( float3 dp, float b ){
    const float r2   = dot(dp,dp);
    const float Rc   = 5.f/b;
    if(r2>(Rc*Rc)){ return (float4){0.f,0.f,0.f,0.f}; }
    const float r    = sqrt(r2);
    const float y    = 1.f-(b*r/5.f);
    const float dy   = b;
    const float y2   = y*y;
    const float y4   = y2*y2;
    return  (float4){ dp*(y4*dy/r), y4*y  };
}

inline float4 exp_r_lin8( float3 dp, float b ){
    const float r2   = dot(dp,dp);
    const float Rc   = 9.f/b;
    if(r2>(Rc*Rc)){ return (float4){0.f,0.f,0.f,0.f}; }
    const float r    = sqrt(r2);
    const float y    = 1.f-(b*r/9.f);
    const float dy   = b;
    const float y2   = y*y;
    const float y4   = y2*y2;
    const float y8   = y4*y4;
    return  (float4){ dp*(y8*dy/r), y8*y  };
}

inline float4 exp_r_lin16( float3 dp, float b ){
    const float r2   = dot(dp,dp);
    const float Rc   = 17.f/b;
    if(r2>(Rc*Rc)){ return (float4){0.f,0.f,0.f,0.f}; }
    const float r    = sqrt(r2);
    const float y    = 1.f-(b*r/17.f);
    const float dy   = b;
    const float y2   = y*y;
    const float y4   = y2*y2;
    const float y8   = y4*y4;
    const float y16  = y8*y8;
    return  (float4){ dp*(y16*dy/r), y16*y  };
}

inline float4 exp_r_cub4( float3 dp, float4 cpoly, float Rc ){
    const float r2   = dot(dp,dp);
    if(r2>(Rc*Rc)){ return (float4){0.f,0.f,0.f,0.f}; }
    const float r    = sqrt(r2);
    const float y    = cpoly.x + r*( cpoly.y + r*(    cpoly.z + r*    cpoly.w)) ;
    const float dy   =               cpoly.y + r*(2.f*cpoly.z + r*3.f*cpoly.w)  ;
    const float y2   = y*y;
    const float y4   = y2*y2;
    return  (float4){ dp*(y4*-5.f*dy/r), y4*y  };
}

inline float4 getMorse( float3 dp,  float R0, float E0, float b ){
    float r = length( dp );
    float e = exp ( -b*(r-R0));
    float E = E0*        e*(e - 2.f); // Energy
    float F = E0*  2.f*b*e*(e - 1.f); // Force
    return  (float4){ dp*(F/r), E };
}

/**
 * @brief Morse potential using a linear approximation of the exponential term (n=5).
 * @param dp float3 distance vector between particles.
 * @param R0 float  Equilibrium distance for the potential minimum.
 * @param E0 float  Potential well depth (magnitude of energy at the minimum).
 * @param b  float  Exponential decay constant, controls the width of the well.
 * @return float4   Force vector in xyz, Energy in w.
 */
inline float4 getMorse_lin5( float3 dp, float R0, float E0, float b ){
    const float n = 5.f;
    const float r2 = dot(dp,dp);
    const float Rc = R0 + n/b;
    if (r2 > Rc*Rc) { return (float4)(0.f); }

    const float r_inv = rsqrt(r2);
    const float r     = r2 * r_inv;
    
    // Base of the power approximation
    const float y = 1.f - (b/n)*(r - R0);
    // The cutoff above handles r > Rc, where y would be negative.

    // p = y^5
    const float y2 = y*y;
    const float y4 = y2*y2;
    const float p  = y4*y;

    // Energy: E = E0 * (p^2 - 2p)
    const float E = E0 * p * (p - 2.f);

    // Force: F = -dE/dr * (dp/r)
    // dp/dr = -b * y^4
    const float dpdr = -b * y4;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;

    return (float4)(dp * f_over_r, E);
}

/**
 * @brief Morse potential using a linear approximation of the exponential term (n=9).
 */
inline float4 getMorse_lin9( float3 dp, float R0, float E0, float b ){
    const float n = 9.f;
    const float r2 = dot(dp,dp);
    const float Rc = R0 + n/b;
    if (r2 > Rc*Rc) { return (float4)(0.f); }

    const float r_inv = rsqrt(r2);
    const float r     = r2 * r_inv;
    
    const float y = 1.f - (b/n)*(r - R0);

    // p = y^9
    const float y2 = y*y;
    const float y4 = y2*y2;
    const float y8 = y4*y4;
    const float p  = y8*y;

    const float E = E0 * p * (p - 2.f);

    const float dpdr = -b * y8;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;

    return (float4)(dp * f_over_r, E);
}

/**
 * @brief Morse potential using a linear approximation of the exponential term (n=17).
 */
inline float4 getMorse_lin17( float3 dp, float R0, float E0, float b ){
    const float n = 17.f;
    const float r2 = dot(dp,dp);
    const float Rc = R0 + n/b;
    if (r2 > Rc*Rc) { return (float4)(0.f); }

    const float r_inv = rsqrt(r2);
    const float r     = r2 * r_inv;
    
    const float y = 1.f - (b/n)*(r - R0);

    // p = y^17
    const float y2  = y*y;
    const float y4  = y2*y2;
    const float y8  = y4*y4;
    const float y16 = y8*y8;
    const float p   = y16*y;

    const float E = E0 * p * (p - 2.f);

    const float dpdr = -b * y16;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;

    return (float4)(dp * f_over_r, E);
}

/**
 * @brief Morse potential using a cubic polynomial approximation of the exponential term (n=5).
 * @param dp    float3 distance vector between particles.
 * @param R0    float  Equilibrium distance for the potential minimum.
 * @param E0    float  Potential well depth.
 * @param cpoly float4 Coefficients (c0, c1, c2, c3) of the polynomial in (r-R0).
 * @param Rc    float  Cutoff radius for the interaction.
 * @return float4   Force vector in xyz, Energy in w.
 */
inline float4 getMorse_cub5( float3 dp, float R0, float E0, float4 cpoly, float Rc ){
    const float r2 = dot(dp,dp);
    if (r2 > Rc*Rc) { return (float4)(0.f); }
    //const float r_inv = rsqrt(r2);
    //const float r     = r2 * r_inv;
    const float r     = sqrt(r2);
    const float r_inv = 1/r;
    const float x     = r - R0;
    // Polynomial y = P(x) and its derivative dy/dx, where x = r-R0
    const float y    = cpoly.x + x*( cpoly.y + x*(    cpoly.z + x*    cpoly.w ));
    const float dydr =               cpoly.y + x*(2.f*cpoly.z + x*3.f*cpoly.w);
    // p = y^5
    const float y2 = y*y;
    const float y4 = y2*y2;
    const float p  = y4*y;
    // Energy: E = E0 * (p^2 - 2p)
    const float E = E0 * p * (p - 2.f);
    // Force: F = -dE/dr * (dp/r)
    // dp/dr = 5 * y^4 * dy/dr
    const float dpdr     =  5.f * y4 * dydr;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;
    return (float4)(dp * f_over_r, E);
}

/**
 * @brief Morse potential using a cubic polynomial approximation (n=9).
 */
inline float4 getMorse_cub9( float3 dp, float R0, float E0, float4 cpoly, float Rc ){
    const float r2 = dot(dp,dp);
    if (r2 > Rc*Rc) { return (float4)(0.f); }

    const float r_inv = rsqrt(r2);
    const float r     = r2 * r_inv;
    const float x     = r - R0;

    const float y    = cpoly.x + x*( cpoly.y + x*( cpoly.z + x*cpoly.w ));
    const float dydr =             cpoly.y + x*(2.f*cpoly.z + x*3.f*cpoly.w);

    // p = y^9
    const float y2 = y*y;
    const float y4 = y2*y2;
    const float y8 = y4*y4;
    const float p  = y8*y;
    
    const float E = E0 * p * (p - 2.f);

    // dp/dr = 9 * y^8 * dy/dr
    const float dpdr = 9.f * y8 * dydr;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;

    return (float4)(dp * f_over_r, E);
}

/**
 * @brief Morse potential using a cubic polynomial approximation (n=17).
 */
inline float4 getMorse_cub17( float3 dp, float R0, float E0, float4 cpoly, float Rc ){
    const float r2 = dot(dp,dp);
    if (r2 > Rc*Rc) { return (float4)(0.f); }

    const float r_inv = rsqrt(r2);
    const float r     = r2 * r_inv;
    const float x     = r - R0;

    const float y    = cpoly.x + x*( cpoly.y + x*( cpoly.z + x*cpoly.w ));
    const float dydr =             cpoly.y + x*(2.f*cpoly.z + x*3.f*cpoly.w);

    // p = y^17
    const float y2 = y*y;
    const float y4 = y2*y2;
    const float y8 = y4*y4;
    const float y16= y8*y8;
    const float p  = y16*y;
    
    const float E = E0 * p * (p - 2.f);

    // dp/dr = 17 * y^16 * dy/dr
    const float dpdr = 17.f * y16 * dydr;
    const float f_scalar = -2.f * E0 * (p - 1.f) * dpdr;
    const float f_over_r = f_scalar * r_inv;

    return (float4)(dp * f_over_r, E);
}




// evaluate non-covalent interaction force and energy for Lennard-Jones (Q) and Coulomb interactions of charges (Q) and hydrogen bond correction (pseudo-charges H), damping R2damp is used to avoid singularity at r=0
inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
    // ---- Electrostatic (damped Coulomb potential)
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 +  R2damp);              // inverse distance squared and damped
    float   Ec    =  COULOMB_CONST*REQ.z*sqrt( ir2_ ); // Ec = Q1*Q2/sqrt(r^2+R2damp)
    // --- Lennard-Jones and Hydrogen bond correction
    float  ir2 = 1.f/r2;          // inverse distance squared
    float  u2  = REQ.x*REQ.x*ir2; // u2 = (R0/r)^2
    float  u6  = u2*u2*u2;        // u6 = (R0/r)^6
    float vdW  = u6*REQ.y;        // vdW = E0*(R0/r)^6
    float E    =       (u6-2.f)*vdW     + Ec  ;     // E = E0*(R0/r)^6 - E0*(R0/r)^12 + Q1*Q2/sqrt(r^2+R2damp)
    float fr   = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;  // fr = -12*E0*( (R0/r)^8/r + 12*E0*(R0/r)^14) - Q1*Q2/(r^2+R2damp)^1.5
    return  (float4){ dp*fr, E };
}

inline float4 getMorseQH( float3 dp,  float4 REQH, float K, float R2damp ){
    float r2    = dot(dp,dp);
    float ir2_  = 1/(r2+R2damp);
    float r     = sqrt( r2   );
    float ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    float e     = exp ( K*(r-REQH.x));
    //double e2    = e*e;
    //double fMors =  E0*  2*K*( e2 -   e ); // Morse
    //double EMors =  E0*      ( e2 - 2*e );
    float   Ae  = REQH.y*e;
    float fMors = Ae*  2*K*(e - 1); // Morse
    float EMors = Ae*      (e - 2);
    float Eel   = COULOMB_CONST*REQH.z*ir_;
    float fr    = fMors/r - Eel*ir2_ ;
    return  (float4){ dp*fr, EMors+Eel };
}

// evaluate damped Coulomb potential and force 
inline float4 getCoulomb( float3 dp, float R2damp ){
    // ---- Electrostatic
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 + R2damp);
    float   E    = COULOMB_CONST*sqrt( ir2_ );
    return  (float4){ dp*-E*ir2_, E };
}

// limit force magnitude to fmax
float3 limnitForce( float3 f, float fmax ){
    float fr2 = dot(f,f);                         // force magnitude squared
    if( fr2>(fmax*fmax) ){ f*=(fmax/sqrt(fr2)); } // if force magnitude is larger than fmax we scale it down to fmax
    return f;
}

float4 getR4repulsion( float3 d, float R, float Rcut, float A ){
    // we use R4blob(r) = A * (1-r^2)^2
    // such that at distance r=R we have force f = fmax
    // f = -dR4blob/dr = 4*A*r*(1-r^2) = fmax
    // A = fmax/(4*R*(1-R^2))
    float R2    = R*R;
    float R2cut = Rcut*Rcut;
    float r2 = dot(d,d);
    if( r2>R2cut ){ 
        return (float4){0.0f,0.0f,0.0f,0.0f};
    }else if( r2>R2 ){ 
        float mr2 = R2cut-r2;
        float fr = A*mr2;
        return (float4){ d*(-4*fr), fr*mr2 };
    }else{
        float mr2 = R2cut-R2;
        float fr  = A*mr2;
        float r    = sqrt(r2);
        float fmax = 4*R*fr;
        return (float4){ d* (-fmax/r), fmax*(R-r) + fr*mr2 };
    }
}


// ======================================================================
// ======================================================================
//                           MMFF kernells
// ======================================================================
// ======================================================================


// ======================================================================
//                     cleanForceMMFFf4()
// ======================================================================
// Clean forces on atoms and neighbors to prepare for next forcefield evaluation
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void cleanForceMMFFf4(
    const int4        nDOFs,       // 1
    __global float4*  aforce,      // 5
    __global float4*  fneigh       // 6
){
    const int natoms=nDOFs.x;
    const int nnode =nDOFs.y;
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    const int nvec = natoms+nnode;

    const int iav = iG + iS*nvec;
    const int ian = iG + iS*nnode;

    aforce[iav]=float4Zero;
    //aforce[iav]=(float4){iG,iS,iav,0.0};

    //if(iav==0){ printf("OCL::cleanForceMMFFf4() iS %i nG %i nS %i \n", iS, nG, nS );}
    //if(iG==0){ for(int i=0;i<(natoms+nnode);i++ ){printf("cleanForceMMFFf4[%i](%g,%g,%g)\n",i,aforce[i].x,aforce[i].y,aforce[i].z);} }
    if(iG<nnode){ 
        const int i4 = ian*4;
        fneigh[i4+0]=float4Zero;
        fneigh[i4+1]=float4Zero;
        fneigh[i4+2]=float4Zero;
        fneigh[i4+3]=float4Zero;
    }
    //if(iG==0){ printf( "OCL::updateAtomsMMFFf4() END\n" ); }
}

// ======================================================================
//                          getMMFFf4()
// ======================================================================

// 1.  getMMFFf4() - computes bonding interactions between atoms and nodes and its neighbors (max. 4 neighbors allowed), the resulting forces on atoms are stored "aforce" array and recoil forces on neighbors are stored in "fneigh" array
//                   kernel run over all atoms and all systems in parallel to exploit GPU parallelism
//__attribute__((reqd_work_group_size(1,1,1)))
//void func_getMMFFf4(
__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    __global float4*  REQs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,
    const int npbc,
    const int bSubtractVdW
){

    const int iG = get_global_id (0);   // intex of atom   (iG<nAtoms)
    const int iS = get_global_id (1);   // index of system (iS<nS)
    //const int nG = get_global_size(0);
    //const int nS = get_global_size(1);  // number of systems
    //const int iL = get_local_id  (0);
    //const int nL = get_local_size(0);
    const int nAtoms=nDOFs.x;  // number of atoms in the system
    const int nnode =nDOFs.y;  // number of nodes in the system
    //const int nvec  = nAtoms+nnode;

    if(iG>=nnode) return;

    const int i0a   = iS*nAtoms;         // index of first atom      in the system
    const int i0n   = iS*nnode;          // index of first node atom in the system
    const int i0v   = iS*(nAtoms+nnode); // index of first vector    in the system ( either atom or pi-orbital )
    
    const int iaa = iG + i0a;  // index of current atom (either node or capping atom)
    const int ian = iG + i0n;  // index of current node atom
    const int iav = iG + i0v;  // index of current vector ( either atom or pi-orbital )

    #define NNEIGH 4
    
    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds (h.xyz) and inverse bond lengths (h.w)
    float3  fbs[4];              // force on neighbor sigma    (fbs[i] is sigma recoil force on i-th neighbor)
    float3  fps[4];              // force on neighbor pi       (fps[i] is pi    recoil force on i-th neighbor)
    float3  fa  = float3Zero;    // force on center atom positon 

    float E=0;                   // Total Energy of this atom
    // ---- Params
    const int4   ng  = neighs[iaa];    // neighboring atoms
    const float3 pa  = apos[iav].xyz;  // position of current atom
    const float4 par = apars[ian];     // (xy=s0_ss,z=ssK,w=piC0 ) forcefield parameters for current atom

    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getMMFFf4(): natoms %3i nnode %3i nSys %3i \n", nAtoms, nnode, get_global_size(1) ); }
    // if((iG==iGdbg)&&(iS==iSdbg)){
    //     printf("OCL::getMMFFf4(): natoms %3i nnode %3i \n", nAtoms, nnode );
    //     for(int ia=0; ia<nnode; ia++){
    //         int4   ng=neighs[iaa+ia];
    //         float4 pi=apos[iav+ia];
    //         float4 bk=bLs[ian];
    //         float4 bK=bKs[ian];
    //         float4 Ks=Ksp[ian];
    //         float4 Kp=Kpp[ian];
    //         float4 apar=apars[ian];
    //         //printf("OCL::getMMFFf4(): ia %3i: pos=(%10.5f,%10.5f,%10.5f) ngs=(%3i,%3i,%3i,%3i) bLs=(%10.5f,%10.5f,%10.5f,%10.5f) bKs=(%10.5f,%10.5f,%10.5f,%10.5f) Ks=(%10.5f,%10.5f,%10.5f,%10.5f) Kp=(%10.5f,%10.5f,%10.5f,%10.5f) apar=(%10.5f,%10.5f,%10.5f,%10.5f)\n", 
    //         //   ia, pi.x, pi.y, pi.z, ng.x, ng.y, ng.z, ng.w, bk.x, bk.y, bk.z, bk.w, bK.x, bK.y, bK.z, bK.w, Ks.x, Ks.y, Ks.z, Ks.w, Kp.x, Kp.y, Kp.z, Kp.w, apar.x, apar.y, apar.z, apar.w );
    //         printf("OCL::getMMFFf4(): ia %3i: pos=(%10.5f,%10.5f,%10.5f) ngs=(%3i,%3i,%3i,%3i) bLs=(%10.5f,%10.5f,%10.5f,%10.5f) bKs=(%10.5f,%10.5f,%10.5f,%10.5f) apar=(%10.5f,%10.5f,%10.5f,%10.5f)\n", 
    //            ia, pi.x, pi.y, pi.z, ng.x, ng.y, ng.z, ng.w, bk.x, bk.y, bk.z, bk.w, bK.x, bK.y, bK.z, bK.w, apar.x, apar.y, apar.z, apar.w );
    //     } 
    // }

    // Temp Arrays
    const int*   ings  = (int*  )&ng; // neighboring atoms, we cast it to int[] to be index it in for loop


    const float   ssC0   = par.x*par.x - par.y*par.y;                      // cos(2) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf , where ang0 is equilibrium angle
    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }   // clear recoil forces on neighbors

    float3 f1,f2;         // working forces 

    { // ========= BONDS - here we evaluate pairwise interactions of node atoms with its 4 neighbors 

        float3  fpi = float3Zero;                // force on pi-orbital
        const int4   ngC = neighCell[iaa];       // neighboring atom cell index
        const float3 hpi = apos[iav+nAtoms].xyz; // direction of pi-orbital
        const float4 vbL = bLs[ian];             // bond lengths
        const float4 vbK = bKs[ian];             // bond stiffness
        const float4 vKs = Ksp[ian];             // stiffness of sigma-pi othogonalization 
        const float4 vKp = Kpp[ian];             // stiffness of pi-pi    alignment

        const int*   ingC  = (int*  )&ngC;   // neighboring atom cell index (we cast it to int[] to be index it in for loop)
        const float* bL    = (float*)&vbL;   // bond lengths
        const float* bK    = (float*)&vbK;   // bond stiffness
        const float* Kspi  = (float*)&vKs;   // stiffness of sigma-pi othogonalization
        const float* Kppi  = (float*)&vKp;   // stiffness of pi-pi    alignment

        const int ipbc0 = iS*npbc;  // index of first PBC shift for current system

        for(int i=0; i<NNEIGH; i++){  // loop over 4 neighbors
            float4 h;                 // direction vector of bond
            const int ing  = ings[i]; // index of i-th neighbor node atom
            const int ingv = ing+i0v; // index of i-th neighbor vector
            const int inga = ing+i0a; // index of i-th neighbor atom
            if(ing<0) break;
            
            // --- Compute bond direction vector and inverse bond length
            h.xyz    = apos[ingv].xyz - pa;  // direction vector of bond
            { // shift bond to the proper PBC cell
                int ic  = ingC[i];                  // index of i-th neighbor cell
                h.xyz  += pbc_shifts[ipbc0+ic].xyz; // shift bond to the proper PBC cell
            }
            float  l = length(h.xyz);  // compute bond length
            h.w      = 1.f/l;           // store ivnerse bond length
            h.xyz   *= h.w;            // normalize bond direction vector
            hs[i]    = h;              // store bond direction vector and inverse bond length

            //float epp = 0; // pi-pi    energy
            //float esp = 0; // pi-sigma energy

            // --- Evaluate bond-length stretching energy and forces
            if(iG<ing){  
                float elb = evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1; E+=elb;  // harmonic bond stretching, fa is force on center atom, fbs[i] is recoil force on i-th neighbor, 
                    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getMMFFf4(): bond-length:    iG %3i ing %3i elb %10.5f f(%10.5f,%10.5f,%10.5f) l %10.5f bL %10.5f bK %10.5f \n", iG, ing, elb, f1.x, f1.y, f1.z, l, bL[i], bK[i] ); }

                // pi-pi alignment interaction            
                float kpp = Kppi[i];
                if( (ing<nnode) && (kpp>1.e-6f) ){   // Only node atoms have pi-pi alignemnt interaction
                    float epp = evalPiAling( hpi, apos[ingv+nAtoms].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2; E+=epp;    //   pi-alignment(konjugation), fpi is force on pi-orbital, fps[i] is recoil force on i-th neighbor's pi-orbital
                    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getMMFFf4(): cos(pi,pi):     iG %3i ing %3i epp %10.5f f(%10.5f,%10.5f,%10.5f) c %10.5f kpp %10.5f \n", iG, ing, epp, f1.x, f1.y, f1.z, dot(hpi.xyz,apos[ingv+nAtoms].xyz), kpp ); }
                }
            } 
            
            // pi-sigma othogonalization interaction
            float ksp = Kspi[i];
            if(ksp>1.e-6f){  
                float esp = evalAngCos( (float4){hpi,1.f}, h, ksp, par.w, &f1, &f2 );   fpi+=f1; fa-=f2;  fbs[i]+=f2; E+=esp;    //   pi-planarization (orthogonality), fpi is force on pi-orbital, fbs[i] is recoil force on i-th neighbor   
                //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::etMMFFf4(): cos(pi,sigma):      iG %3i ing %3i esp %10.5f f1(%10.5f,%10.5f,%10.5f) c %10.5f ksp %10.5f par.w %10.5f \n", iG, ing, esp, f1.x, f1.y, f1.z, dot(hpi.xyz,h.xyz), ksp, par.w ); }
            }
        }

        // --- Store Pi-forces                      we store pi-forces here because we don't use them in the angular force evaluation
        const int i4p=(iG + iS*nnode*2 )*4 + nnode*4; // index of first pi-force for current atom
        for(int i=0; i<NNEIGH; i++){
            fneigh[i4p+i] = (float4){fps[i],0}; // store recoil pi-force on i-th neighbor
        }
        aforce[iav+nAtoms]  = (float4){fpi,0};  // store pi-force on pi-orbital of current atom

    }
    
    { //  ============== Angles   - here we evaluate angular interactions between pair of sigma-bonds of node atoms with its 4 neighbors

        for(int i=0; i<NNEIGH; i++){ // loop over first bond
            int ing = ings[i];
            if(ing<0) break;         // if there is no i-th neighbor we break the loop
            const float4 hi = hs[i];
            const int ingv = ing+i0v;
            const int inga = ing+i0a;
            for(int j=i+1; j<NNEIGH; j++){ // loop over second bond
                int jng  = ings[j];
                if(jng<0) break;           // if there is no j-th neighbor we break the loop
                const int jngv = jng+i0v;
                const int jnga = jng+i0a;
                const float4 hj = hs[j];  
                      
                float ea = evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );    // evaluate angular force and energy using cos(angle/2) formulation        
                fa  -= f1+f2;
                E   += ea;
                //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getMMFFf4(): angle():        iG %3i ing %3i jng %3i ea=%10.5f f1(%10.5f,%10.5f,%10.5f) f2(%10.5f,%10.5f,%10.5f) cos %10.5f apar(%10.5f,%10.5f,%10.5f) \n", iG, ing, jng, ea, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z, dot(hi.xyz,hj.xyz), par.x, par.y, par.z ); }

                //if(bSubtractVdW)
                { // Remove non-bonded interactions from atoms that are bonded to common neighbor
                    float4 REQi=REQs[inga];   // non-bonding parameters of i-th neighbor
                    float4 REQj=REQs[jnga];   // non-bonding parameters of j-th neighbor
                    // combine non-bonding parameters of i-th and j-th neighbors using mixing rules
                    float4 REQij;             
                    REQij.x  = REQi.x  + REQj.x;
                    REQij.yz = REQi.yz * REQj.yz; 
                    
                    float3 dp = (hj.xyz/hj.w) - (hi.xyz/hi.w);   // recover vector between i-th and j-th neighbors using stored vectos and inverse bond lengths, this should be faster than dp=apos[jngv].xyz-apos[ingv].xyz; from global memory
                    float4 fij = getLJQH( dp, REQij, 1.0f );     // compute non-bonded interaction between i-th and j-th neighbors using LJQH potential
                    f1 -=  fij.xyz;
                    f2 +=  fij.xyz;
                }

                fbs[i]+= f1;
                fbs[j]+= f2;
            }
        }

    }
    
    // ========= Save results - store forces on atoms and recoil on its neighbors  (pi-forces are already done)
    const int i4 =(iG + iS*nnode*2 )*4;
    //const int i4p=i4+nnode*4;
    for(int i=0; i<NNEIGH; i++){
        fneigh[i4 +i] = (float4){fbs[i],0};
        //fneigh[i4p+i] = (float4){fps[i],0};
    }
    //aforce[iav     ] = (float4){fa ,0}; // If we do  run it as first forcefield
    aforce[iav       ] += (float4){fa ,0};  // If we not run it as first forcefield
    //aforce[iav+nAtoms]  = (float4){fpi,0}; 
    
}

// __kernel void getMMFFf4(
//     const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
//     // Dynamical
//     __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
//     __global float4*  aforce,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
//     __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
//     // parameters
//     __global int4*    neighs,       // 5  [nnode]  neighboring atoms
//     __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
//     __global float4*  REQs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
//     __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
//     __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
//     __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
//     __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
//     __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
//     __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
//     __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
//     __global float4*  pbc_shifts,
//     const int npbc,
//     const int bSubtractVdW
// ){
//     func_getMMFFf4( nDOFs, apos, aforce, fneigh, neighs, neighCell, REQs, apars, bLs, bKs, Ksp, Kpp, lvecs, ilvecs, pbc_shifts, npbc, bSubtractVdW );
// }

// ======================================================================
//                           getNonBond()
// ======================================================================

// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
//void func_getNonBond(
__kernel void getNonBond(
    const int4        nDOFs,        // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,         // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams
    //,     // 9 // Grid-Force-Field parameters
    //__local float4*   LATOMS,
    //__local float4*   LCLJS
){
    // we use local memory to store atomic position and parameters to speed up calculation, the size of local buffers should be equal to local workgroup size
    //__local float4 LATOMS[2];
    //__local float4 LCLJS [2];
    //__local float4 LATOMS[4];
    //__local float4 LCLJS [4];
    //__local float4 LATOMS[8];
    //__local float4 LCLJS [8];
    //__local float4 LATOMS[16];
    //__local float4 LCLJS [16];
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters
    //__local float4 LATOMS[64];
    //__local float4 LCLJS [64];
    //__local float4 LATOMS[128];
    //__local float4 LCLJS [128];

    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=nDOFs.x;  // number of atoms
    const int nnode =nDOFs.y;  // number of node atoms
    //const int nAtomCeil =ns.w;
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)
    //const int i0n = iS*nnode; 
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    //const int ian = iG + i0n;
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array
    
    //const int iS_DBG = 0;
    //const int iG_DBG = 0;

    // NOTE: if(iG>=natoms) we are reading from invalid adress => last few processors produce crap, but that is not a problem
    //       importaint is that we do not write this crap to invalid address, so we put   if(iG<natoms){forces[iav]+=fe;} at the end
    //       we may also put these if(iG<natoms){ .. } around more things, but that will unnecessarily slow down other processors
    //       we need these processors with (iG>=natoms) to read remaining atoms to the local memory.

    //if(iG<natoms){
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;  // PBC is used if any of the PBC dimensions is >0
    //const bool bPBC=false;

    const int4   ng     = neighs   [iaa];  // neighbors indices
    const int4   ngC    = neighCell[iaa];  // neighbors cell indices
    const float4 REQKi  = REQs     [iaa];  // non-bonded parameters
    const float3 posi   = apos     [iav].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe           = float4Zero;  // force on atom

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors for this system
    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getNonBond() getNonBond(): natoms=%i, nnode=%i nSys=%i nPBC=(%i,%i,%i)\n", natoms, nnode, get_global_size(1), nPBC.x, nPBC.y, nPBC.z); }
    // if((iG==iGdbg)&&(iS==iSdbg)){ 
    //     printf("OCL::getNonBond() getNonBond(): natoms=%i, nnode=%i nPBC=(%i,%i,%i)\n", natoms, nnode, nPBC.x, nPBC.y, nPBC.z);
    //     printf("OCL::getNonBond(): lvec.a=(%g,%g,%g) lvec.b=(%g,%g,%g) lvec.c=(%g,%g,%g)\n", lvec.a.x, lvec.a.y, lvec.a.z, lvec.b.x, lvec.b.y, lvec.b.z, lvec.c.x, lvec.c.y, lvec.c.z);
    //     printf("OCL::getNonBond(): GFFParams=(%g,%g,%g,%g) \n", GFFParams.x, GFFParams.y, GFFParams.z, GFFParams.w);
    //     for(int i=0; i<natoms; i++){
    //         float4 pi = apos[i];
    //         int4   ng = neighs[i];
    //         int4   ngC = neighCell[i];
    //         float4 REQKi = REQs[i];
    //         printf("OCL::getNonBond() atom %i: ng=(%i,%i,%i,%i), ngC=(%i,%i,%i,%i), REQKi=(%10.5f,%10.5f,%10.5f|%10.5f), posi=(%10.5f,%10.5f,%10.5f,%10.5f)\n", i, ng.x, ng.y, ng.z, ng.w, ngC.x, ngC.y, ngC.z, ngC.w, REQKi.x, REQKi.y, REQKi.z, REQKi.w, pi.x, pi.y, pi.z, pi.w);
    //     }   
    // }

    //if(iG==0){ printf("GPU[iS=%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f} \n", iS, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );  }

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQs[i].x,REQs[i].y,REQs[i].z ); }

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;   // shift of PBC image 0
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                      // shift of PBC image in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                      // shift of PBC image in the outer loop
    //}
    /*
    if((iG==iG_DBG)&&(iS==iS_DBG)){ 
        printf( "OCL::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) bPBC=%i nPBC(%i,%i,%i)\n", natoms,nnode,nvec, nS,nG,nL, bPBC, nPBC.x,nPBC.y,nPBC.z ); 
        for(int i=0; i<natoms; i++){
            printf( "GPU a[%i] ", i);
            printf( "p{%6.3f,%6.3f,%6.3f} ", atoms[i0v+i].x,atoms[i0v+i].y,atoms[i0v+i].z  );
            printf( "ng{%i,%i,%i,%i} ", neighs[i0a+i].x,neighs[i0a+i].y,neighs[i0a+i].z,neighs[i0a+i].w );
            printf( "ngC{%i,%i,%i,%i} ", neighCell[i0a+i].x,neighCell[i0a+i].y,neighCell[i0a+i].z,neighCell[i0a+i].w );
            printf( "\n");
        }
    }
    */

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){     // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;             // index of atom in local memory
        if(i<natoms){                  // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = apos [i+i0v]; // read atom position to local memory 
            LCLJS [iL] = REQs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                //if( (j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //if( (ji==1)&&(iG==0) )printf( "2 non-bond[%i,%i] bBonded %i\n",iG,ji,bBonded );

                if(bPBC){         // ===== if PBC is used, we need to loop over all PBC images of the atom
                    int ipbc=0;   // index of PBC image
                    dp += shift0; // shift to PBC image 0
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "GPU ipbc %i(%i,%i) shift(%7.3g,%7.3g,%7.3g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z ); 
                            // Without these IF conditions if(bBonded) time of evaluation reduced from 61 [ms] to 51[ms]
                            if( !( bBonded &&(                     // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                      ((ja==ng.x)&&(ipbc==ngC.x))    // check if this PBC image is not the same as one of the bonded atoms
                                    ||((ja==ng.y)&&(ipbc==ngC.y))  // i.e. if ja is neighbor of iG, and ipbc is its neighbor cell index then we skip this interaction
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "ipbc %i(%i,%i) shift(%g,%g,%g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z ); 
                                float4 fij = getLJQH( dp, REQK, R2damp );  // calculate non-bonded force between atoms using LJQH potential
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                fe += fij;
                            }
                            ipbc++; 
                            dp    += lvec.a.xyz; 
                        }
                        dp    += shift_a;
                    }
                }else                                                  // ===== if PBC is not used, it is much simpler
                if( !bBonded ){  
                    float4 fij = getLJQH( dp, REQK, R2damp ); 
                    fe += fij;
                    //if((iG==iGdbg)&&(iS==iSdbg)){   printf("OCL::getNonBond(): ia,ja %3i %3i aj(%10.5f,%10.5f,%10.5f) dp( %10.5f | %10.5f,%10.5f,%10.5f)  fij( %10.5f,%10.5f,%10.5f|%10.5f)\n", iG, ja, aj.x,aj.y,aj.z, length(dp), dp.x, dp.y, dp.z, fij.x, fij.y, fij.z, fij.w); }
                }  // if atoms are not bonded, we calculate non-bonded interaction between them
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    if(iG<natoms){
        //if(iS==0){ printf( "OCL::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        aforce[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)
        //aforce[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //aforce[iav] = fe*(-1.f);
    }
}



__kernel void scanNonBond(
    const    int     n,       // 1  number of points
    const    float4  REQH,    // 2  non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global float4* pos,   // 3  [n]positions of points
    __global float4* force, // 4  [n]forces on points
    const    float8  ffpar    // 5  parameters specific to the potential function used
){
    for(int i=0; i<n; i++){
        float4 fij;
        const float3 dp = pos[i].xyz;
        //fij=getForce(dp,REQH,ffpar.x);
        //<<<GET_FORCE_NONBOND   // this line will be replaced python pre-processor
        force[i] = fij;
    }
}

#define WG_scanNonBond2 32
__kernel void scanNonBond2(
    const int         n,     // 1  number of points
    const float4      REQH0, // 2  non-bonded parameters of test atom (RvdW,EvdW,QvdW,Hbond)
    __global float4*  pos,   // 3  [n] positions of points
    __global float4*  force, // 4  [n] forces on points
    const int         na,    // 5  number of atoms
    __global float4*  apos,  // 6  [na] postions of atoms
    __global float4*  REQs,  // 7  [na] non-bonded parameters of atoms (RvdW,EvdW,QvdW,Hbond)
    const float8      ffpar  // 8  parameters specific to the potential function used
){
    __local float4 lPos[WG_scanNonBond2];          // cached atom positions
    __local float4 lPar[WG_scanNonBond2];          // cached atom parameters
    const int    lid = get_local_id(0);  // 0 … WG-1
    const int    gid = get_global_id(0); // 0 … n-1
    const float3 p   = pos[gid].xyz;     // position of the test point
    float4       f   = (float4)(0.0f);
    // tile atoms through shared memory
    for(int il0=0; il0<na; il0+=WG_scanNonBond2){
        int ia=il0+lid;      // global atom index
        if(ia<na){
            lPos[lid]=apos[ia];
            lPar[lid]=REQs[ia];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // loop over the current tile
        #pragma unroll           // optional
        for (int j=0; j<WG_scanNonBond2; ++j) {
            int ja=il0+j;
            if(ia<na){
                const float3 dp  = lPos[j].xyz-p;
                float4 REQH      = lPar[j];
                REQH.x  +=REQH0.x;
                REQH.yzw*=REQH0.yzw;
                float4 fij;
                //fij=getForce(dp,REQH,ffpar.x);
                //<<<GET_FORCE_NONBOND   // this line will be replaced python pre-processor
                f+=fij;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(gid<n){
        force[gid]=f;
    }
}

__kernel void scanNonBond2PBC( 
    const int         n,      // 1  number of points
    const float4      REQH0,  // 2  non-bonded parameters of test atom (RvdW,EvdW,QvdW,Hbond)
    __global float4*  pos,    // 3  [n] positions of points
    __global float4*  force,  // 4  [n] forces on points
    const int         na,     // 5  number of atoms
    __global float4*  apos,   // 6  [na] postions of atoms
    __global float4*  REQs,   // 7  [na] non-bonded parameters of atoms (RvdW,EvdW,QvdW,Hbond)
    const float8      ffpar,  // 8  parameters specific to the potential function used
    const cl_Mat3     lvec,   // 9  lattice vectors for each system
    const int4        nPBC    // 10 number of PBC images in each direction (x,y,z)
){
    __local float4 lPos[WG_scanNonBond2];          // cached atom positions
    __local float4 lPar[WG_scanNonBond2];          // cached atom parameters
    const int    lid = get_local_id(0);  // 0 … WG-1
    const int    gid = get_global_id(0); // 0 … n-1
    const float3 p   = pos[gid].xyz;     // position of the test point
    float4       f   = (float4)(0.0f);
    // tile atoms through shared memory
    for(int il0=0; il0<na; il0+=WG_scanNonBond2){
        int ia=il0+lid;      // global atom index
        if(ia<na){
            lPos[lid]=apos[ia];
            lPar[lid]=REQs[ia];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // loop over the current tile
        #pragma unroll           // optional
        for (int j=0; j<WG_scanNonBond2; ++j) {
            int ja=il0+j;
            if(ia<na){
                const float3 dp0  = lPos[j].xyz-p;
                float4 REQH      = lPar[j];
                REQH.x  +=REQH0.x;
                REQH.yzw*=REQH0.yzw;
                for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                    const float3 dp0x = dp0 + lvec.b.xyz*ix;
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        const float3 dp0y = dp0x + lvec.c.xyz*iy;
                        for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                            const float3 dp = dp0y + lvec.a.xyz*iz;
                            float4 fij;
                            //fij=getForce(dp,REQH,ffpar.x);
                            //<<<GET_FORCE_NONBOND   // this line will be replaced python pre-processor
                            f+=fij;
                        }
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(gid<n){
        force[gid]=f;
    }
}

__kernel void scanNonBond2PBC_2( 
    const int         n,      // 1  number of points
    const float4      REQH0,  // 2  non-bonded parameters of test atom (RvdW,EvdW,QvdW,Hbond)
    __global float4*  pos,    // 3  [n] positions of points
    __global float4*  force,  // 4  [n] forces on points
    const int         na,     // 5  number of atoms
    __global float4*  apos,   // 6  [na] postions of atoms
    __global float4*  REQs,   // 7  [na] non-bonded parameters of atoms (RvdW,EvdW,QvdW,Hbond)
    const float8      ffpar,  // 8  parameters specific to the potential function used
    const cl_Mat3     lvec,   // 9  lattice vectors for each system
    const int4        nPBC    // 10 number of PBC images in each direction (x,y,z)
){
    __local float4 lPos[WG_scanNonBond2];          // cached atom positions
    __local float4 lPar[WG_scanNonBond2];          // cached atom parameters
    const int    lid = get_local_id(0);  // 0 … WG-1
    const int    gid = get_global_id(0); // 0 … n-1
    const float3 p   = pos[gid].xyz;     // position of the test point
    float4       f   = (float4)(0.0f);
    // tile atoms through shared memory
    for(int il0=0; il0<na; il0+=WG_scanNonBond2){
        int ia=il0+lid;      // global atom index
        if(ia<na){
            lPos[lid]=apos[ia];
            lPar[lid]=REQs[ia];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // loop over the current tile
        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
            const float3 p0x = p + lvec.b.xyz*ix;
            for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                const float3 p0y = p0x + lvec.c.xyz*iy;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    const float3 p0 = p0y + lvec.a.xyz*iz;
                    #pragma unroll           // optional
                    for (int j=0; j<WG_scanNonBond2; ++j) {
                        int ja=il0+j;
                        if(ia<na){
                            const float3 dp  = lPos[j].xyz-p0;
                            float4 REQH      = lPar[j];
                            REQH.x  +=REQH0.x;
                            REQH.yzw*=REQH0.yzw;
                            float4 fij;
                            //fij=getForce(dp,REQH,ffpar.x);
                            //<<<GET_FORCE_NONBOND   // this line will be replaced python pre-processor
                            f+=fij;
                        }
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(gid<n){
        force[gid]=f;
    }
}



__kernel void getNonBond_template(
    const int4        nDOFs,        // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,         // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams
    //,     // 9 // Grid-Force-Field parameters
    //__local float4*   LATOMS,
    //__local float4*   LCLJS
){
    // we use local memory to store atomic position and parameters to speed up calculation, the size of local buffers should be equal to local workgroup size
    //__local float4 LATOMS[2];
    //__local float4 LCLJS [2];
    //__local float4 LATOMS[4];
    //__local float4 LCLJS [4];
    //__local float4 LATOMS[8];
    //__local float4 LCLJS [8];
    //__local float4 LATOMS[16];
    //__local float4 LCLJS [16];
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters
    //__local float4 LATOMS[64];
    //__local float4 LCLJS [64];
    //__local float4 LATOMS[128];
    //__local float4 LCLJS [128];

    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=nDOFs.x;  // number of atoms
    const int nnode =nDOFs.y;  // number of node atoms
    //const int nAtomCeil =ns.w;
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)
    //const int i0n = iS*nnode; 
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    //const int ian = iG + i0n;
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array
    
    //const int iS_DBG = 0;
    //const int iG_DBG = 0;

    // NOTE: if(iG>=natoms) we are reading from invalid adress => last few processors produce crap, but that is not a problem
    //       importaint is that we do not write this crap to invalid address, so we put   if(iG<natoms){forces[iav]+=fe;} at the end
    //       we may also put these if(iG<natoms){ .. } around more things, but that will unnecessarily slow down other processors
    //       we need these processors with (iG>=natoms) to read remaining atoms to the local memory.

    //if(iG<natoms){
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;  // PBC is used if any of the PBC dimensions is >0
    //const bool bPBC=false;

    const int4   ng     = neighs   [iaa];  // neighbors indices
    const int4   ngC    = neighCell[iaa];  // neighbors cell indices
    const float4 REQKi  = REQs     [iaa];  // non-bonded parameters
    const float3 posi   = apos     [iav].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe           = float4Zero;  // force on atom

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors for this system
    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::getNonBond() getNonBond(): natoms=%i, nnode=%i nSys=%i nPBC=(%i,%i,%i)\n", natoms, nnode, get_global_size(1), nPBC.x, nPBC.y, nPBC.z); }
    // if((iG==iGdbg)&&(iS==iSdbg)){ 
    //     printf("OCL::getNonBond() getNonBond(): natoms=%i, nnode=%i nPBC=(%i,%i,%i)\n", natoms, nnode, nPBC.x, nPBC.y, nPBC.z);
    //     printf("OCL::getNonBond(): lvec.a=(%g,%g,%g) lvec.b=(%g,%g,%g) lvec.c=(%g,%g,%g)\n", lvec.a.x, lvec.a.y, lvec.a.z, lvec.b.x, lvec.b.y, lvec.b.z, lvec.c.x, lvec.c.y, lvec.c.z);
    //     printf("OCL::getNonBond(): GFFParams=(%g,%g,%g,%g) \n", GFFParams.x, GFFParams.y, GFFParams.z, GFFParams.w);
    //     for(int i=0; i<natoms; i++){
    //         float4 pi = apos[i];
    //         int4   ng = neighs[i];
    //         int4   ngC = neighCell[i];
    //         float4 REQKi = REQs[i];
    //         printf("OCL::getNonBond() atom %i: ng=(%i,%i,%i,%i), ngC=(%i,%i,%i,%i), REQKi=(%10.5f,%10.5f,%10.5f|%10.5f), posi=(%10.5f,%10.5f,%10.5f,%10.5f)\n", i, ng.x, ng.y, ng.z, ng.w, ngC.x, ngC.y, ngC.z, ngC.w, REQKi.x, REQKi.y, REQKi.z, REQKi.w, pi.x, pi.y, pi.z, pi.w);
    //     }   
    // }

    //if(iG==0){ printf("GPU[iS=%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f} \n", iS, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );  }

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQs[i].x,REQs[i].y,REQs[i].z ); }

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;   // shift of PBC image 0
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                      // shift of PBC image in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                      // shift of PBC image in the outer loop
    //}
    /*
    if((iG==iG_DBG)&&(iS==iS_DBG)){ 
        printf( "OCL::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) bPBC=%i nPBC(%i,%i,%i)\n", natoms,nnode,nvec, nS,nG,nL, bPBC, nPBC.x,nPBC.y,nPBC.z ); 
        for(int i=0; i<natoms; i++){
            printf( "GPU a[%i] ", i);
            printf( "p{%6.3f,%6.3f,%6.3f} ", atoms[i0v+i].x,atoms[i0v+i].y,atoms[i0v+i].z  );
            printf( "ng{%i,%i,%i,%i} ", neighs[i0a+i].x,neighs[i0a+i].y,neighs[i0a+i].z,neighs[i0a+i].w );
            printf( "ngC{%i,%i,%i,%i} ", neighCell[i0a+i].x,neighCell[i0a+i].y,neighCell[i0a+i].z,neighCell[i0a+i].w );
            printf( "\n");
        }
    }
    */

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){     // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;             // index of atom in local memory
        if(i<natoms){                  // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = apos [i+i0v]; // read atom position to local memory 
            LCLJS [iL] = REQs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                //if( (j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //if( (ji==1)&&(iG==0) )printf( "2 non-bond[%i,%i] bBonded %i\n",iG,ji,bBonded );

                if(bPBC){         // ===== if PBC is used, we need to loop over all PBC images of the atom
                    int ipbc=0;   // index of PBC image
                    dp += shift0; // shift to PBC image 0
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "GPU ipbc %i(%i,%i) shift(%7.3g,%7.3g,%7.3g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z ); 
                            // Without these IF conditions if(bBonded) time of evaluation reduced from 61 [ms] to 51[ms]
                            if( !( bBonded &&(                     // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                      ((ja==ng.x)&&(ipbc==ngC.x))    // check if this PBC image is not the same as one of the bonded atoms
                                    ||((ja==ng.y)&&(ipbc==ngC.y))  // i.e. if ja is neighbor of iG, and ipbc is its neighbor cell index then we skip this interaction
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                // float4 fij = getLJQH( dp, REQK, R2damp );
                                float4 fij;
                                //<<<---GET_FORCE_NONBOND   // this line will be replaced python pre-processor
                                fe += fij;
                            }
                            ipbc++; 
                            dp    += lvec.a.xyz; 
                        }
                        dp    += shift_a;
                    }
                }else   // ===== if PBC is not used, it is much simpler
                if( !bBonded ){  
                    //float4 fij = getLJQH( dp, REQK, R2damp ); 
                    float4 fij;
                    //<<<---GET_FORCE_NONBOND   // this line will be replaced python pre-processor
                    fe += fij;
                    //if((iG==iGdbg)&&(iS==iSdbg)){   printf("OCL::getNonBond(): ia,ja %3i %3i aj(%10.5f,%10.5f,%10.5f) dp( %10.5f | %10.5f,%10.5f,%10.5f)  fij( %10.5f,%10.5f,%10.5f|%10.5f)\n", iG, ja, aj.x,aj.y,aj.z, length(dp), dp.x, dp.y, dp.z, fij.x, fij.y, fij.z, fij.w); }
                }  // if atoms are not bonded, we calculate non-bonded interaction between them
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    if(iG<natoms){
        //if(iS==0){ printf( "OCL::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        aforce[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)
        //aforce[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //aforce[iav] = fe*(-1.f);
    }
}



// __kernel void getNonBond(
//     const int4        nDOFs,        // 1 // (natoms,nnode) dimensions of the system
//     // Dynamical
//     __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
//     __global float4*  aforce,       // 3 // forces on atoms
//     // Parameters
//     __global float4*  REQs,         // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
//     __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
//     __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
//     __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
//     const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
//     const float4      GFFParams     // 9 // Grid-Force-Field parameters
// ){
//     __local float4 LATOMS[32];
//     __local float4 LCLJS[32];
//     func_getNonBond( nDOFs, apos, aforce, REQs, neighs, neighCell, lvecs, nPBC, GFFParams, LATOMS, LCLJS );
// }


// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

// unsigned int hash_wang(unsigned int bits) {
//     //unsigned int bits = __float_as_int(value);
//     bits = (bits ^ 61) ^ (bits >> 16);
//     bits *= 9;
//     bits = bits ^ (bits >> 4);
//     bits *= 0x27d4eb2d;
//     bits = bits ^ (bits >> 15);
//     return bits;
// }

// float hashf_wang( float val, float xmin, float xmax) {
//     //return ( (float)(bits)*(2147483647.0f );
//     return (((float)( hash_wang(  __float_as_int(val) ) )) * 4.6566129e-10 )  *(xmax-xmin)+ xmin;
// }

// Assemble recoil forces from neighbors and  update atoms positions and velocities 
//__attribute__((reqd_work_group_size(1,1,1)))
//  void func_updateAtomsMMFFf4(
__kernel void updateAtomsMMFFf4(    
    const int4        nDOFs,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms 
    __global float4*  aforce,       // 4 // forces on atoms
    __global float4*  cvf,          // 5 // damping coefficients for velocity and force
    __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
    __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
    __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
){
    const int natoms=nDOFs.x;           // number of atoms
    const int nnode =nDOFs.y;           // number of node atoms
    const int nMaxSysNeighs = nDOFs.z;  // max number of inter-system interactions; if <0 shwitch inter system interactions off
    const int nvec  = natoms+nnode; // number of vectors (atoms+node atoms)
    const int iG = get_global_id  (0); // index of atom

    if(iG>=nvec) return;

    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0); // number of atoms
    const int nS = get_global_size(1); // number of systems

    //const int ian = iG + iS*nnode; 
    const int iaa = iG + iS*natoms;  // index of atom in atoms array
    const int iav = iG + iS*nvec;    // index of atom in vectors array

    const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit)
    const float4 TDrive = TDrives[iS];

    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("OCL::updateAtomsMMFFf4() natoms %i nvec %i nSys %i MDpars(%g,%g,%g,%g) \n", natoms, nvec, get_global_size(1), MDpars.x,MDpars.y,MDpars.z,MDpars.w); }
    // if((iG==iGdbg)&&(iS==iSdbg)){ 
    //     printf("OCL updateAtomsMMFFf4() natoms %i nvec %i MDpars(%g,%g,%g,%g) \n", natoms, nvec, MDpars.x,MDpars.y,MDpars.z,MDpars.w);  
    //     for(int is=0; is<nS; is++){
    //         //printf( "OCL::TDrives[%i](%g,%g,%g,%g)\n", i, TDrives[i].x,TDrives[i].y,TDrives[i].z,TDrives[i].w );
    //         //printf( "OCL::bboxes[%i](%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", is, bboxes[is].a.x,bboxes[is].a.y,bboxes[is].a.z,   bboxes[is].b.x,bboxes[is].b.y,bboxes[is].b.z,   bboxes[is].c.x,bboxes[is].c.y,bboxes[is].c.z );
    //         for(int ia=0; ia<natoms; ia++){
    //             int ic = ia+is*natoms;
    //             if(constr[ia+is*natoms].w>0) printf( "OCL  sys[%i]atom[%i] constr(%g,%g,%g|K=%g) constrK(%g,%g,%g|%g)\n", is, ia, constr[ic].x,constr[ic].y,constr[ic].z,constr[ic].w,   constrK[ic].x,constrK[ic].y,constrK[ic].z,constrK[ic].w  );
    //         }
    //     }
    // }

    //if((iG==iGdbg)&&(iS==iSdbg))printf( "updateAtomsMMFFf4() natoms=%i nnode=%i nvec=%i nG %i iS %i/%i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS, MDpars.x, MDpars.y, MDpars.z );
    // if((iG==iGdbg)&&(iS==iSdbg)){
    //     int i0a = iS*natoms;
    //     for(int i=0; i<natoms; i++){
    //         printf( "GPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   );
    //     }
    // }
    
    if(iG>=(natoms+nnode)) return; // make sure we are not out of bounds of current system

    //aforce[iav] = float4Zero;

    float4 fe      = aforce[iav]; // force on atom or pi-orbital
    const bool bPi = iG>=natoms;  // is it pi-orbital ?
    
    // ------ Gather Forces from back-neighbors

    int4 ngs = bkNeighs[ iav ]; // back neighbors indices

    //if(iS==5)printf( "iG,iS %i %i ngs %i,%i,%i,%i \n", iG, iS, ngs.x,ngs.y,ngs.z,ngs.w );
    //if( (iS==0)&&(iG==0) ){ printf( "GPU:fe.1[iS=%i,iG=%i](%g,%g,%g,%g) \n", fe.x,fe.y,fe.z,fe.w ); }

    // sum all recoil forces from back neighbors   - WARRNING : bkNeighs must be properly shifted on CPU by adding offset of system iS*nvec*4
    if(ngs.x>=0){ fe += fneigh[ngs.x]; } // if neighbor index is negative it means that there is no neighbor, so we skip it
    if(ngs.y>=0){ fe += fneigh[ngs.y]; }
    if(ngs.z>=0){ fe += fneigh[ngs.z]; }
    if(ngs.w>=0){ fe += fneigh[ngs.w]; }

    // ---- Limit Forces - WARRNING : Github_Copilot says: this is not the best way to limit forces, because it can lead to drift, better is to limit forces in the first forcefield run (best is NBFF) 
    float Flimit = 10.0f;
    float fr2 = dot(fe.xyz,fe.xyz);  // squared force
    if( fr2 > (Flimit*Flimit) ){  fe.xyz*=(Flimit/sqrt(fr2)); }  // if force is too big, we scale it down to Flimit

    // =============== FORCE DONE
    //aforce[iav] = fe;           // store force before limit
    //aforce[iav] = float4Zero;   // clean force   : This can be done in the first forcefield run (best is NBFF)

    //if((iG==iGdbg)&&(iS==iSdbg)){ printf( "OCL updateAtomsMMFFf4() fe[iS=%3i,iG=%3i](%16.8f,%16.8f,%16.8f|%16.8f) \n", iS,iG, fe.x,fe.y,fe.z,fe.w ); }

    
    // =============== DYNAMICS

    float4 ve = avel[iav]; // velocity of atom or pi-orbital
    float4 pe = apos[iav]; // position of atom or pi-orbital

    // -------- Fixed Atoms and Bounding Box
    if(iG<natoms){                  // only atoms have constraints, not pi-orbitals
        // ------- bboxes
        const cl_Mat3 B = bboxes[iS];
        // if(B.c.x>0.0f){ if(pe.x<B.a.x){ fe.x+=(B.a.x-pe.x)*B.c.x; }else if(pe.x>B.b.x){ fe.x+=(B.b.x-pe.x)*B.c.x; }; }
        // if(B.c.y>0.0f){ if(pe.y<B.a.y){ fe.y+=(B.a.y-pe.y)*B.c.y; }else if(pe.y>B.b.y){ fe.y+=(B.b.y-pe.y)*B.c.y; }; }
        if(B.c.z>0.0f){ if(pe.z<B.a.z){ fe.z+=(B.a.z-pe.z)*B.c.z; }else if(pe.z>B.b.z){ fe.z+=(B.b.z-pe.z)*B.c.z; }; }
        // ------- constrains
        float4 cons = constr[ iaa ]; // constraints (x,y,z,K)
        if( cons.w>0.f ){            // if stiffness is positive, we have constraint
            float4 cK = constrK[ iaa ];
            cK = max( cK, (float4){0.0f,0.0f,0.0f,0.0f} );
            const float3 fc = (cons.xyz - pe.xyz)*cK.xyz;
            fe.xyz += fc; // add constraint force
            //if(iS==0){printf( "OCL::constr[ia=%i|iS=%i] (%g,%g,%g|K=%g) fc(%g,%g,%g) cK(%g,%g,%g)\n", iG, iS, cons.x,cons.y,cons.z,cons.w, fc.x,fc.y,fc.z , cK.x, cK.y, cK.z ); }
        }
    }

    // -------- Inter system interactions
    if( nMaxSysNeighs>0 ){
        for(int i=0; i<nMaxSysNeighs; i++){
            const int j     = iS*nMaxSysNeighs + i;
            const int    jS = sysneighs[j];
            const float4 bj = sysbonds [j];
            const float4 pj = apos[jS*nvec + iG];
            float3 d        = pj.xyz - pe.xyz;
            float  l = length( d );
            if      (l<bj.x){
                d*=(l-bj.x)*bj.z/l;  // f = dx*kPress
            }else if(l>bj.y){
                d*=(bj.y-l)*bj.w/l;  // f = dx*kTens
            }
            fe.xyz += d;
        }
    }
    
    const bool bDrive = TDrive.y > 0.0f;

    // ------ Move (Leap-Frog)
    if(bPi){ // if pi-orbital, we need to make sure that it has unit length
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght, 
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }else{
        // Thermal driving  - Langevin thermostat, see C++ MMFFsp3_loc::move_atom_Langevin()
        if( bDrive ){ // if gamma>0
            fe.xyz    += ve.xyz * -TDrive.y ;  // damping,  check the untis  ... cdamp/dt = gamma
            //const float3 rnd = (float3){ hashf_wang(ve.x+TDrive.w,-1.0f,1.0f),hashf_wang(ve.y+TDrive.w,-1.0f,1.0f),hashf_wang(ve.z+TDrive.w,-1.0f,1.0f)};
            __private float3 ix; 
            // + (float3){TDrive.w,TDrive.w,TDrive.w}
            //const float4 rnd = fract( (ve*541547.1547987f + TDrive.wwww), &ix )*2.f - (float4){1.0f,1.0f,1.0f,1.0f};  // changes every frame
            const float3 rvec = (float3){  // random vector depending on the index
                (((iG+136  + (int)(1000.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f, 
                (((iG+778  + (int)(1013.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f,
                (((iG+4578 + (int)( 998.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f
            };
            //const float3 rnd = fract( ( rvec + TDrive.www)*12.4565f, &ix )*2.f - (float3){1.0f,1.0f,1.0f};
            const float3 rnd = sin( ( rvec + TDrive.www )*124.4565f );
            //if(iS==3){  printf( "atom[%i] seed=%g rvec(%g,%g,%g) rnd(%g,%g,%g) \n", iG, TDrive.w, rvec.x,rvec.y,rvec.z, rnd.x,rnd.y,rnd.z ); }
            fe.xyz    += rnd.xyz * sqrt( 2*const_kB*TDrive.x*TDrive.y/MDpars.x );
        }
    }
    cvf[iav] += (float4){ dot(fe.xyz,fe.xyz),dot(ve.xyz,ve.xyz),dot(fe.xyz,ve.xyz), 0.0f };    // accumulate |f|^2 , |v|^2  and  <f|v>  to calculate damping coefficients for FIRE algorithm outside of this kernel
    //if(!bDrive){ ve.xyz *= MDpars.z; } // friction, velocity damping
    ve.xyz *= MDpars.z;             // friction, velocity damping
    ve.xyz += fe.xyz*MDpars.x;      // acceleration
    pe.xyz += ve.xyz*MDpars.x;      // move
    //ve     *= 0.99f;              // friction, velocity damping
    //ve.xyz += fe.xyz*0.1f;        // acceleration
    //pe.xyz += ve.xyz*0.1f;        // move
    if(bPi){        // if pi-orbital, we need to make sure that it has unit length
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;    // This seems to be needed, not sure why ?????
    avel[iav] = ve;   // store velocity
    apos[iav] = pe;   // store position
    

    /*
    //------ Move Gradient-Descent
    if(bPi){ fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz ); } // subtract forces  component which change pi-orbital lenght
    //pe.xyz += fe.xyz*MDpars.x*0.01f;
    //pe.xyz += fe.xyz*MDpars.x*0.01f;
    pe.xyz += fe.xyz*0.02f;
    //if(bPi){ pe.xyz=normalize(pe.xyz); }
    pe.w=0;  // This seems to be needed, not sure why ?????
    apos[iav] = pe;
    */
    
    //if(iG==0){ printf( "OCL::updateAtomsMMFFf4() END\n" ); }
    
}

// __kernel void updateAtomsMMFFf4(
//     const int4        nDOFs,            // 1 // (natoms,nnode) dimensions of the system
//     __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
//     __global float4*  avel,         // 3 // velocities of atoms 
//     __global float4*  aforce,       // 4 // forces on atoms
//     __global float4*  cvf,          // 5 // damping coefficients for velocity and force
//     __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
//     __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
//     __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
//     __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
//     __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
//     __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
//     __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
//     __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
//     __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
// ){
//     func_updateAtomsMMFFf4( nDOFs, apos, avel, aforce, cvf, fneigh, bkNeighs, constr, constrK, MDparams, TDrives, bboxes, sysneighs, sysbonds );
// }

__kernel void runMD(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms 
    __global float4*  aforce,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    __global float4*  REQs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge     

    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,
    const int npbc,
    const int bSubtractVdW,

    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams,     // 9 // Grid-Force-Field parameters

    __global float4*  cvf,          // 5 // damping coefficients for velocity and force
    __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
    __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
){

   const int iG = get_global_id(0);
   const int iS = get_global_id(1);
   //if((iG==iGdbg)&&(iS==iSdbg)){ printf( "OCL runMD() nDOFs(%i,%i,%i,%i) \n", nDOFs.x,nDOFs.y,nDOFs.z,nDOFs.w ); }

    //__local float4 LATOMS[32];
    //__local float4 LCLJS[32];
    for(int istep=0; istep<nDOFs.w; istep++){
        //if((iG==iGdbg)&&(iS==iSdbg)){ printf( "OCL::runMD() --- iter %i nDOFs(%i,%i,%i,%i) \n", istep, nDOFs.x,nDOFs.y,nDOFs.z,nDOFs.w ); }
        //func_getMMFFf4(
        getMMFFf4(
            nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
            apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
            aforce,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
            fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
            neighs,       // 5  [nnode]  neighboring atoms
            neighCell,    // 5  [nnode]  neighboring atom  cell index
            REQs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge 
            apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
            bLs,          // 8  [nnode]  bond length    between node and each neighbor
            bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
            Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
            Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
            lvecs,        // 12 lattice vectors         for each system
            ilvecs,       // 13 inverse lattice vectors for each system
            pbc_shifts,
            npbc,
            bSubtractVdW
        );
        barrier(CLK_GLOBAL_MEM_FENCE);
        //func_getNonBond(
        getNonBond(
            nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
            apos,         // 2  [natoms]     positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
            aforce,      // 3 // forces on atoms
            REQs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
            neighs,      // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
            neighCell,   // 6 // neighbors cell indices ( to ignore interactions between bonded atoms )
            lvecs,        // 7 // lattice vectors for each system
            nPBC,         // 8 // number of PBC images in each direction (x,y,z)
            GFFParams
            //,    // 9 // Grid-Force-Field parameters
            //LATOMS,
            //LCLJS
        );
        barrier(CLK_GLOBAL_MEM_FENCE);
        //func_updateAtomsMMFFf4(
        updateAtomsMMFFf4(
            nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
            apos,         // 2  [natoms]     positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
            avel,         // 3 // velocities of atoms 
            aforce,       // 4 // forces on atoms
            cvf,          // 5 // damping coefficients for velocity and force
            fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
            bkNeighs,     // 7 // back neighbors indices (for recoil forces)
            constr,       // 8 // constraints (x,y,z,K) for each atom
            constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
            MDparams,     // 10 // MD parameters (dt,damp,Flimit)
            TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
            bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
            sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
            sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
        );
        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}





/*

// ======================================================================
//                     printOnGPU()
// ======================================================================
// Print atoms and forces on GPU
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void printOnGPU(
    const int4        nDOFs,        // 1
    const int4        mask,         // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
){
    const int natoms=nDOFs.x;
    const int nnode =nDOFs.y; 
    const int isys  =nDOFs.z; 
    const int nvec  = natoms+nnode;
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    //const int ian = iG + iS*nnode; 
    const int iaa = iG + iS*natoms; 
    const int iav = iG + iS*nvec;
    //const int iG_DBG = 0;
    const int iG_DBG = 1;
    
    printf( "#### OCL::printOnGPU(isys=%i) natoms=%i nnode=%i nG,nS(%i,%i) \n", isys,  natoms,nnode,   nS,nG );
    if(mask.x){
        for(int i=0; i<natoms; i++){
            int ia=i + isys*nvec;
            printf( "GPU[%i=%i] ", i, ia );
            //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[ia].x, bkNeighs[ia].y, bkNeighs[ia].z, bkNeighs[ia].w );
            printf( "aforce{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[ia].x, aforce[ia].y, aforce[ia].z, aforce[ia].w );
            //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[ia].x, avel[ia].y, avel[ia].z, avel[ia].w );
            printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[ia].x, apos[ia].y, apos[ia].z, apos[ia].w );
            printf(  "constr{%6.3f,%6.3f,%6.3f,%6.3f} ", constr[ia].x, constr[ia].y, constr[ia].z, constr[ia].w );
            printf( "\n" );
        }
    }
    if(mask.y){
        for(int i=0; i<nnode; i++){
            int i1=i+natoms + isys*nvec;
            printf( "GPU[%i=%i] ", i, i1 );
            printf(  "fpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i1].x, aforce[i1].y, aforce[i1].z, aforce[i1].w );
            //printf(  "vpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z, avel[i1].w );
            printf(   "pipos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z, apos[i1].w );
            printf( "\n" );
        }
    }
    if(mask.z){
        for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
            int i0 = isys*nvec*8;
            int i1=i0 + i*4+j;
            int i2=i0 + i1 + nnode*4;
            printf( "GPU[%i,%i] ", i, j );
            printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i1].x, fneigh[i1].y, fneigh[i1].z, fneigh[i1].w );
            printf( "fneighpi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i2].x, fneigh[i2].y, fneigh[i2].z, fneigh[i2].w );
            printf( "\n" );
        }}
    }
    
}

*/


// ======================================================================
//         B-spline Interpolation Functions using buffer
// ======================================================================

inline int modulo(const int i, const int m) {
    int result = i % m;
    if (result < 0) {
        result += m;
    }
    return result;
}

inline int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
    // iqs[0] = (int4)(0, 1,   2,   3  );
    // iqs[1] = (int4)(0, 1,   2,   3-n);
    // iqs[2] = (int4)(0, 1,   2-n, 3-n);
    // iqs[3] = (int4)(0, 1-n, 2-n, 3-n);
}

inline int4 choose_inds_pbc(const int i, const int n, const int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
    if(i>=(n-3)){ 
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const int4 d = iqs[ii];
        return (int4){ i+d.x, i+d.y, i+d.z, i+d.w }; 
    }
    return (int4){ i, i+1, i+2, i+3 };
}


inline float4 basis(float u) {
    const float inv6 = 1.0f / 6.0f;
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        inv6 * t * t * t,
        inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f),
        inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f),
        inv6 * u2 * u
    );
}

inline float4 dbasis(float u) {
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        -0.5f * t * t,
        0.5f * (3.0f * u2 - 4.0f * u),
        0.5f * (-3.0f * u2 + 2.0f * u + 1.0f),
        0.5f * u2
    );
}

inline float2 fe1Dcomb(__global const float4* E, const float4 C, const float4 p, const float4 d) {
    const float4 cs = (float4)(dot(C, E[0]), dot(C, E[1]), dot(C, E[2]), dot(C, E[3]));
    return (float2)(dot(p, cs), dot(d, cs));
}

inline float3 fe2d_comb(int nz, __global const float4* E, int4 di, const float4 C, const float4 pz, const float4 dz, const float4 by, const float4 dy) {
    const float2 fe0 = fe1Dcomb(E + di.x, C, pz, dz);
    const float2 fe1 = fe1Dcomb(E + di.y, C, pz, dz);
    const float2 fe2 = fe1Dcomb(E + di.z, C, pz, dz);
    const float2 fe3 = fe1Dcomb(E + di.w, C, pz, dz);
    
    return (float3)(
        fe0.x * dy.x + fe1.x * dy.y + fe2.x * dy.z + fe3.x * dy.w,
        fe0.y * by.x + fe1.y * by.y + fe2.y * by.z + fe3.y * by.w,
        fe0.x * by.x + fe1.x * by.y + fe2.x * by.z + fe3.x * by.w
    );
}

inline float4 fe3d_pbc_comb(const float3 u, const int3 n, __global const float4* Es, const float4 PLQH, __local const int4* xqis, __local int4* yqis) {
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;

    if ((iz < 1) || (iz >= n.z - 2)) {
        return (float4)(0.12356f, 0.0f, 0.0f, 0.0f);
    }

    ix = modulo(ix-1, n.x);
    iy = modulo(iy-1, n.y);

    const int nyz = n.z * n.y;
    // int4 qx = xqis[ix%4] * nyz;
    // int4 qy = yqis[iy%4] * n.z;

    int4 qx = choose_inds_pbc( ix, n.x, xqis );
    //const int4 qx = choose_inds_pbc( ix, n.x, xqis )*nyz;
    const int4 qy = choose_inds_pbc( iy, n.y, yqis )*n.z;

    const float4 bz = basis(tz);
    const float4 dz = dbasis(tz);
    const float4 by = basis(ty);
    const float4 dy = dbasis(ty);
    
    const int i0 = (iz - 1) + n.z * (iy + n.y * ix);

    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) n(%i,%i,%i) \n", u.x,u.y,u.z, ix,iy,iz, n.x,n.y,n.z );
    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) qx(%i,%i,%i,%i) nyz=%i\n", u.x,u.y,u.z, ix,iy,iz, qx.x,qx.y,qx.z,qx.w, nyz );
    qx*=nyz;
    
    //return (float4){ 0.0f, 0.0f, 0.0f, dot(PLQH, Es[ i0 ])  };
    //return  Es[i0] ;

    float3 E1 = fe2d_comb(n.z, Es + (i0 + qx.x), qy, PLQH, bz, dz, by, dy);
    float3 E2 = fe2d_comb(n.z, Es + (i0 + qx.y), qy, PLQH, bz, dz, by, dy);
    float3 E3 = fe2d_comb(n.z, Es + (i0 + qx.z), qy, PLQH, bz, dz, by, dy);
    float3 E4 = fe2d_comb(n.z, Es + (i0 + qx.w), qy, PLQH, bz, dz, by, dy);
    const float4 bx = basis(tx);
    const float4 dx = dbasis(tx);
    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))
    );
    
}

int pos2int( const float3 pos, const int4 grid_ns ){
    return (int)(pos.z) + grid_ns.z * ( (int)(pos.y) + grid_ns.y * (int)(pos.x) );
}

__kernel void sampleGrid(
    const int4 nDOFs,                      // 1
    __global float4*  apos,                // 2
    __global float4*  aforce,              // 3
    __global float4*  REQs,                // 4
    const float4  GFFParams,               // 5
    __global float4* BsplinePLQH,          // 6
    const float4   grid_invStep,           // 7
    const int4     grid_ns,                // 8
    const float4   grid_p0,                 // 9
    const float4   MDparams,               // 10
    int      nstep                   // 11
){
    const int iG = get_global_id  (0);
    const int iL = get_local_id   (0);
    const int nG = get_global_size(0);
    const int np = nDOFs.x;

    //float3 dz = (float3){ 0.0f, 0.0f, 0.1f };

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    // const float4 REQ        = REQs [iG];
    // const float3 posi       = apos[iG].xyz;
    // const float  R2damp     = GFFParams.x*GFFParams.x;
    // const float  alphaMorse = GFFParams.y;

    //aforce[iG] = fe3d_pbc_comb(posi, grid_ns.xyz, BsplinePLQH, (float4){0.0f,0.0f,0.0f,1.0f});
    //aforce[iG] = BsplinePLQH[ pos2int(posi.xyz,grid_ns) ];





    __local int4 xqs[4];
    __local int4 yqs[4];

    if      (iL<4){             xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
    else if (iL<8){ int i=iL-4; yqs[i ]=make_inds_pbc(grid_ns.y,i ); };
    //const float3 inv_dg = 1.0f / grid_d.xyz;
    barrier(CLK_LOCAL_MEM_FENCE);

    const int iG_DBG = 0;

    //if(iG==iG_DBG){   printf( "GPU::sampleGrid() grid_ns(%i,%i,%i) grid_invStep(%g,%g,%g) grid_p0(%g,%g,%g) \n", grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x, grid_invStep.y, grid_invStep.z, grid_p0.x, grid_p0.y, grid_p0.z ); }
    //if(iG==iG_DBG){   printf( "GPU::sampleGrid() xqs(%i,%i,%i,%i) \n", xqs[0],xqs[1],xqs[2],xqs[3] ); yqs(%i,%i,%i,%i) \n", xqs[0],xqs[1],xqs[2],xqs[3], yqs[0],yqs[1],yqs[2],yqs[3] ); }
    // if(iG==iG_DBG){   
    //     printf( "GPU::sampleGrid() xqs[0](%i,%i,%i,%i) \n", xqs[0].x,xqs[0].y,xqs[0].z,xqs[0].w ); 
    //     printf( "GPU::sampleGrid() xqs[1](%i,%i,%i,%i) \n", xqs[1].x,xqs[1].y,xqs[1].z,xqs[1].w ); 
    //     printf( "GPU::sampleGrid() xqs[2](%i,%i,%i,%i) \n", xqs[2].x,xqs[2].y,xqs[2].z,xqs[2].w ); 
    //     printf( "GPU::sampleGrid() xqs[3](%i,%i,%i,%i) \n", xqs[3].x,xqs[3].y,xqs[3].z,xqs[3].w ); 
    //     printf( "GPU::sampleGrid() yqs[0](%i,%i,%i,%i) \n", yqs[0].x,yqs[0].y,yqs[0].z,yqs[0].w );
    //     printf( "GPU::sampleGrid() yqs[1](%i,%i,%i,%i) \n", yqs[1].x,yqs[1].y,yqs[1].z,yqs[1].w );
    //     printf( "GPU::sampleGrid() yqs[2](%i,%i,%i,%i) \n", yqs[2].x,yqs[2].y,yqs[2].z,yqs[2].w );
    //     printf( "GPU::sampleGrid() yqs[3](%i,%i,%i,%i) \n", yqs[3].x,yqs[3].y,yqs[3].z,yqs[3].w );
    // }


    //const float ej = exp( GFFParams.y * REQ.x ); // exp(-alphaMorse*RvdW) pre-factor for factorized Morse potential
    //const float4 PLQH = (float4){ ej*ej*REQ.y, ej*REQ.y,  REQ.z, 0.0f };
    const float4 PLQH = (float4){ 1.0f, 1.0f, 1.0f, 1.0f };

    const float dt    =      MDparams.x;
    const float cdamp = 1.0f-MDparams.y;
    //nstep = 1;

    //if(iG==iG_DBG){   printf( "GPU::sampleGrid() np=%i nstep=%i dt=%g cdamp=%g grid_ns(%i,%i,%i,%i) \n", np, nstep, dt, cdamp, grid_ns.x,grid_ns.y,grid_ns.z,grid_ns.w ); }

    float3 posi  = apos[iG].xyz;
    float4 fe    = (float4){0.0f,0.0f,0.0f,0.0f};
    float3 vel   = (float3){0.0f,0.0f,0.0f};

    // float3 u = (posi - grid_p0.xyz) * grid_invStep.xyz;
    // fe      = fe3d_pbc_comb(u, grid_ns.xyz, BsplinePLQH, PLQH, xqs, yqs);
    // fe.xyz *= -grid_invStep.xyz;

    for(int i=0; i<nstep; i++){
        float3 u  = (posi - grid_p0.xyz) * grid_invStep.xyz;
        //u        +=(float3){0.5f,0.5f,0.5f};
        fe  = fe3d_pbc_comb(u, grid_ns.xyz, BsplinePLQH, PLQH, xqs, yqs);
        vel       *= cdamp;
        vel       += fe.xyz*dt;
        posi      += vel   *dt;
    }

    aforce[iG]    = fe;

    //float4 fg = BsplinePLQH[ pos2int(u,grid_ns) ];

    // if((iG==iG_DBG)&&(iS==iS_DBG)){   
    //     printf( "GPU[%i] grid_ns(%i,%i,%i) grid_invStep(%g,%g,%g) grid_p0(%g,%g,%g) \n", iG, grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x,grid_invStep.y,grid_invStep.z, grid_p0.x,grid_p0.y,grid_p0.z ); 
    //     //printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8e,%16.8e,%16.8e,%16.8e) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w );   
    //     //float4 feg = read_imagef(BsplinePLQH_tex, sampler_bspline, (int4)(u.x, -5,     u.z, 0));
    //     //printf("GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) feg(%16.8e,%16.8e,%16.8e,%16.8e)\n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, feg.x,feg.y,feg.z,feg.w);
    // }

    //fes[iG] = fe;

    // if(iG==0){ 
    //     printf( "GPU::sampleGridFF() np=%i R2damp=%g aMorse=%g grid_ns(%i,%i,%i,%i) \n", np, R2damp, alphaMorse, grid_ns.x,grid_ns.y,grid_ns.z,grid_ns.w ); 
    //     for(int i=0; i<np; i++){
    //         float4 pi = apos[i];
    //         int idx = pos2int(pi.xyz,grid_ns);
    //         float4 fe = BsplinePLQH[idx];
    //         printf( "GPU::sampleGridFF() %i idx=%8i p(%16.8f,%16.8f,%16.8f)  force(%16.8f,%16.8f,%16.8f|%16.8f) \n", i, idx, pi.x,pi.y,pi.z, fe.x,fe.y,fe.z,fe.w ); 
    //     }; 
    // }
// NOTE: https://registry.khronos.org/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
// CLK_ADDRESS_REPEAT - out-of-range image coordinates are wrapped to the valid range. This address mode can only be used with normalized coordinates. If normalized coordinates are not used, this addressing mode may generate image coordinates that are undefined.
}


// ======================================================================
//                           getNonBond_GridFF_Bspline()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system, and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by texture interpolation
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
// NOTE: this modified version of getNonBond() just by including Grid-Force-Field (GFF) forces
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF_Bspline(
    const int4 nDOFs,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  apos,        // 2 // positions of atoms
    __global float4*  aforce,      // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW,EvdW,Q,H)
    // GridFF
    __global float4*  BsplinePLQH,   // 10 // Grid-Force-Field (GFF) for Pauli repulsion
    const int4     grid_ns,         // 11 // origin of the grid
    const float4   grid_invStep,    // 12 // origin of the grid
    const float4   grid_p0          // 13 // origin of the grid
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms 
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=nDOFs.x;         // number of atoms in the system
    const int nnode =nDOFs.y;         // number of nodes in the system
    const int nvec  =natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system
    
    const float4 REQKi = REQs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = apos    [iav].xyz;       // position of the atom
    float4 fe          = float4Zero;              // forces on the atom

    const int iS_DBG = 0;
    const int iG_DBG = 0;

    // =================== Non-Bonded interaction ( molecule-molecule )

    { // insulate nbff

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "GPU::getNonBond_GridFF_Bspline() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){ 
    //     //printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); 
    //     for(int i=0; i<natoms; i++){
    //         printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) neighs(%4i,%4i,%4i,%4i) REQ(%16.8f,%16.8f,%16.8f,%16.8f,%16.8f) \n", i, apos[iav+i].x,apos[iav+i].y,apos[iav+i].z, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQs[i].x,REQs[i].y,REQs[i].z,REQs[i].w );
    //     }
    // }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom

    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = apos [i+i0v];      // load positions  of atoms into local memory
        LCLJS [iL] = REQs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory 
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )   
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK.x  +=REQKi.x;              // mixing of RvdW radii
                REQK.yz *=REQKi.yz;             // mixing of EvdW and Q
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom            
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){      
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ? 
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    } // insulate nbff

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory 

    // ========== Molecule-Grid interaction with GridFF using tricubic Bspline ================== (see. kernel sample3D_comb() in GridFF.cl

    __local int4 xqs[4];
    __local int4 yqs[4];
    { // insulate gridff
        if      (iL<4){             xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
        else if (iL<8){ int i=iL-4; yqs[i ]=make_inds_pbc(grid_ns.y,i ); };
        //const float3 inv_dg = 1.0f / grid_d.xyz;
        barrier(CLK_LOCAL_MEM_FENCE);

        const float ej = exp( GFFParams.y * REQKi.x ); // exp(-alphaMorse*RvdW) pre-factor for factorized Morse potential
        const float4 PLQH = (float4){ 
            ej*ej*REQKi.y,                   // prefactor London dispersion (attractive part of Morse potential)
            ej*   REQKi.y,                   // prefactor Pauli repulsion   (repulsive part of Morse potential)
            REQKi.z,
            0.0f
        };
        //const float3 p = ps[iG].xyz;
        const float3 u = (posi - grid_p0.xyz) * grid_invStep.xyz;

        float4 fg = fe3d_pbc_comb(u, grid_ns.xyz, BsplinePLQH, PLQH, xqs, yqs);        
        fg.xyz *= -grid_invStep.xyz;
        //fe += fg;

        //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() fg(%g,%g,%g|%g) u(%g,%g,%g) posi(%g,%g,%g) grid_invStep(%g,%g,%g)\n", fg.x,fg.y,fg.z,fg.w,  u.x,u.y,u.z, posi.x,posi.y,posi.z, grid_invStep.x, grid_invStep.y, grid_invStep.z  ); }

        //float4 feg = BsplinePLQH[ (int)( u.z + grid_ns.z * (u.y + grid_ns.y * u.x ) ) ];
        //fe += feg;

        if((iG==iG_DBG)&&(iS==iS_DBG)){   
            printf( "GPU[%i] grid_ns(%i,%i,%i) grid_invStep(%g,%g,%g) grid_p0(%g,%g,%g) \n", iG, grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x,grid_invStep.y,grid_invStep.z, grid_p0.x,grid_p0.y,grid_p0.z ); 
            //printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8e,%16.8e,%16.8e,%16.8e) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w );   
            //float4 feg = read_imagef(BsplinePLQH_tex, sampler_bspline, (int4)(u.x, -5,     u.z, 0));
            //printf("GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) feg(%16.8e,%16.8e,%16.8e,%16.8e)\n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, feg.x,feg.y,feg.z,feg.w);
        }

        //fes[iG] = fe;
    }  // insulate gridff


    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) force(%16.8f,%16.8f,%16.8f,%16.8f) BsplinePLQ=%li \n", iG, posi.x,posi.y,posi.z, fe.x,fe.y,fe.z,fe.w, BsplinePLQ );}

    aforce[iav] = fe;        // If we do    run it as first forcefield, in this case we do not need to clear forces before running this forcefield
    //aforce[iav] += fe;     // If we don't run it as first forcefield, we need to add forces to the forces calculated by previous forcefields
    //aforce[iav] = fe*(-1.f);
    

}








// ======================================================================
//         B-spline Interpolation Functions using 3D Texture
// ======================================================================

__constant sampler_t sampler_bspline = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

__constant sampler_t sampler_norm_rep_near = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
//__constant sampler_t sampler_norm_rep_near = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
//__constant sampler_t sampler_norm_rep_near = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

// 1D B-spline interpolation along Z, combining 4 potential components Reads 4 float4 values from texture at (ix, iy, qz[0..3]),  Combines R,G,B,A channels using C_coeffs (PLQH),  Returns (Energy, dEnergy/duz)
inline float2 fe1Dcomb_tex_(
    __read_only image3d_t img,  // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const float4 C,             // Coefficients for combining P,L,Q,H
    const float4 pz,            // B-spline basis functions for z (B0, B1, B2, B3)
    const float4 dz,            // B-spline basis functions for z (B'0, B'1, B'2, B'3) derivatives
    const int3 i      
){
    // Combine Pauli, London, Coulomb, H-bond components using C = (P,L,Q,H) for each grid point
    const float4 cs = (float4)(
        dot(C, read_imagef(img, sampler_bspline, (int4)(i.x, i.y, i.z  ,0) )), 
        dot(C, read_imagef(img, sampler_bspline, (int4)(i.x, i.y, i.z+1,0) )), 
        dot(C, read_imagef(img, sampler_bspline, (int4)(i.x, i.y, i.z+2,0) )), 
        dot(C, read_imagef(img, sampler_bspline, (int4)(i.x, i.y, i.z+3,0) ))
    );
    // Interpolate energy and its derivative w.r.t. u_z using B-spline basis
    return (float2)(
        dot(pz, cs), // Energy
        dot(dz, cs)  // dEnergy/du_z
    );
}


// 2D B-spline interpolation in YZ plane, for a given X-coordinate slice
// Calls fe1Dcomb_tex four times.
// Returns (dEnergy/duy, dEnergy/duz, Energy)
inline float3 fe2d_comb_tex_(
    __read_only image3d_t img,  // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const int  ix,              // Single X integer coordinate for this slice
    const int4 qy,              // 4 integer Y coordinates
    const int  iz,              // pixel index along z-axis
    const float4 C,             // Coefficients for combining P,L,Q,H
    const float4 pz,            // B-spline polynomials for Z
    const float4 dz,            // B-spline polynomials for Z derivative
    const float4 py,            // B-spline polynomials for Y 
    const float4 dy             // B-spline polynomials for Y derivative
){
    // Interpolate along Z for 4 different Y lines (at the given ix)
    const float2 fe0 = fe1Dcomb_tex_(img, C, pz, dz, (int3){ix, qy.x, iz});
    const float2 fe1 = fe1Dcomb_tex_(img, C, pz, dz, (int3){ix, qy.y, iz});
    const float2 fe2 = fe1Dcomb_tex_(img, C, pz, dz, (int3){ix, qy.z, iz});
    const float2 fe3 = fe1Dcomb_tex_(img, C, pz, dz, (int3){ix, qy.w, iz});

    // feN.x is Energy(yN,      uz_interp)
    // feN.y is dEnergy/duz(yN, uz_interp)

    // Interpolate along Y using results from 1D Z-interpolation
    return (float3)(
        dot(dy, (float4)(fe0.x, fe1.x, fe2.x, fe3.x)),     // dEnergy/du_y = sum_j (B'_j(u_y) * Energy(y_j, u_z_interp))
        dot(py, (float4)(fe0.y, fe1.y, fe2.y, fe3.y)),     // dEnergy/du_z = sum_j (B_j(u_y) * dEnergy/du_z(y_j, u_z_interp))
        dot(py, (float4)(fe0.x, fe1.x, fe2.x, fe3.x))      // Energy       = sum_j (B_j(u_y) * Energy(y_j, u_z_interp))
    );
}

// 3D B-spline interpolation for force and energy
// u: normalized coordinates (fractional cell coordinates)
// n: dimensions of the B-spline grid (texture dimensions)
// img: 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
// PLQH: coefficients to combine the 4 potential components
// xqis, yqis: precomputed PBC index patterns for X and Y dimensions
// Returns (dEnergy/dux, dEnergy/duy, dEnergy/duz, Energy)
inline float4 fe3d_pbc_comb_tex_(
    const float3 u,            // normalized coordinates (fractional voxel coordinates)
    const int3 n,              // dimensions of the B-spline grid (texture dimensions)
    __read_only image3d_t img, // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const float4 PLQH,         // coefficients to combine the 4 potential components
    __local const int4* xqis,  // Patterns from make_inds_pbc for x-dim
    __local const int4* yqis   // Patterns from make_inds_pbc for y-dim
) {
    // Integer part of u (knot index preceding the point)
    // Matches original code's floor logic for ix, iy, iz
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    if (u.z < 0) iz--; // Also apply floor logic to z

    // Fractional part of u (position within the cell defined by knot ix, iy, iz)
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;

    // B-spline interpolation requires 4 knots starting from index (i-1).
    // The indices needed are (i-1, i, i+1, i+2).
    const int ix_ = ix - 1;
    const int iy_ = iy - 1;
    const int iz_ = iz - 1;

    // Boundary condition for Z: if iz_raw_knot is too close to edge, return zero.
    // iz_raw_knot must be in [1, n.z - 3] for full 4-knot support (iz_knot_start must be >=0, iz_knot_start+3 must be < n.z).
    if ((iz < 1) || (iz >= n.z - 2)) {  return (float4)(1.2345f, 0.0f, 0.0f, 0.0f); }

    // Apply PBC for X and Y dimensions to get base knot index for choose_inds_pbc_3
    // The base knot index for choose_inds_pbc_3 should be the starting index *after* modulo,
    // i.e., (ix-1) % n.x.
    const int ix_pbc_base = modulo(ix_, n.x);
    const int iy_pbc_base = modulo(iy_, n.y);

    // Get the 4 absolute integer grid indices for X and Y using PBC logic
    // These indices will be used directly in read_imagef
    const int4 qx = choose_inds_pbc_3(ix_pbc_base, n.x, xqis);
    const int4 qy = choose_inds_pbc_3(iy_pbc_base, n.y, yqis);

    // Calculate B-spline basis functions and their derivatives
    const float4 bz = basis(tz);  const float4 dz = dbasis(tz);
    const float4 by = basis(ty);  const float4 dy = dbasis(ty);
    const float4 bx = basis(tx);  const float4 dx = dbasis(tx);

    //return read_imagef(img, smp, (int4)( qx.x, qy.x, iz, 0 ));
    //return read_imagef(img, smp, (int4)( qx.y, qy.y, iz+1, 0 ));

    //float2 efz = fe1Dcomb_tex(img, PLQH, bz, dz, (int3){ix, iy, iz});
    //return (float4){0.0f,0.0f,efz.y,efz.x};

    // Interpolate along YZ for 4 different X planes
    // E#.x = dE/duy, E#.y = dE/duz, E#.z = E, all at (qx.#, u_y_interp, u_z_interp)
    const float3 E1 = fe2d_comb_tex_(img, qx.x, qy, iz, PLQH, bz, dz, by, dy);
    const float3 E2 = fe2d_comb_tex_(img, qx.y, qy, iz, PLQH, bz, dz, by, dy);
    const float3 E3 = fe2d_comb_tex_(img, qx.z, qy, iz, PLQH, bz, dz, by, dy);
    const float3 E4 = fe2d_comb_tex_(img, qx.w, qy, iz, PLQH, bz, dz, by, dy);

    // Interpolate along X using results from 2D YZ-interpolation
    // Result is (dE/dux, dE/duy, dE/duz, E_total)
    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),     // dEnergy/du_x = sum_i (B'_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),     // dEnergy/du_y = sum_i (B_i(u_x) * dEnergy/du_y(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),     // dEnergy/du_z = sum_i (B_i(u_x) * dEnergy/du_z(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))      // Energy       = sum_i (B_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
    );
    

}









// ======================================================================
//         B-spline Interpolation Functions using 3D Texture  New ( using (float4) texture repeat and normalized coordinates )
// ======================================================================

inline float2 fe1Dcomb_tex(
    __read_only image3d_t img,  // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const float4 C,             // Coefficients for combining P,L,Q,H
    const float4 pz,            // B-spline basis functions for z (B0, B1, B2, B3)
    const float4 dz,            // B-spline basis functions for z (B'0, B'1, B'2, B'3) derivatives
    const float3 uc,            // normalized coordinates (fractional cell coordinates)
    const float duz             // duz = 1.0f/n.z ( shift by one voxel in z-direction)
){
    // Combine Pauli, London, Coulomb, H-bond components using C = (P,L,Q,H) for each grid point
    const float4 cs = (float4)(
        dot(C, read_imagef(img, sampler_norm_rep_near, (float4)(uc.x, uc.y, uc.z          , 0.0f ) )), 
        dot(C, read_imagef(img, sampler_norm_rep_near, (float4)(uc.x, uc.y, uc.z+duz*1.0f , 0.0f ) )), 
        dot(C, read_imagef(img, sampler_norm_rep_near, (float4)(uc.x, uc.y, uc.z+duz*2.0f , 0.0f ) )), 
        dot(C, read_imagef(img, sampler_norm_rep_near, (float4)(uc.x, uc.y, uc.z+duz*3.0f , 0.0f ) ))
    );
    // Interpolate energy and its derivative w.r.t. u_z using B-spline basis
    return (float2)(
        dot(pz, cs), // Energy
        dot(dz, cs)  // dEnergy/du_z
    );
}

inline float3 fe2d_comb_tex(
    __read_only image3d_t img,  // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const float3 uc,            // normalized coordinates (fractional cell coordinates)
    const float4 C,             // Coefficients for combining P,L,Q,H
    const float4 pz,            // B-spline polynomials for Z
    const float4 dz,            // B-spline polynomials for Z derivative
    const float4 py,            // B-spline polynomials for Y 
    const float4 dy,            // B-spline polynomials for Y derivative
    const float3 dc
){
    // Interpolate along Z for 4 different Y lines (at the given ix)
    const float2 fe0 = fe1Dcomb_tex(img, C, pz, dz, (float3){uc.x, uc.y           , uc.z }, dc.z );
    const float2 fe1 = fe1Dcomb_tex(img, C, pz, dz, (float3){uc.x, uc.y+dc.y      , uc.z }, dc.z );
    const float2 fe2 = fe1Dcomb_tex(img, C, pz, dz, (float3){uc.x, uc.y+dc.y*2.0f , uc.z }, dc.z );
    const float2 fe3 = fe1Dcomb_tex(img, C, pz, dz, (float3){uc.x, uc.y+dc.y*3.0f , uc.z }, dc.z );
    // Interpolate along Y using results from 1D Z-interpolation
    return (float3)(
        dot(dy, (float4)(fe0.x, fe1.x, fe2.x, fe3.x)),     // dEnergy/du_y = sum_j (B'_j(u_y) * Energy(y_j, u_z_interp))
        dot(py, (float4)(fe0.y, fe1.y, fe2.y, fe3.y)),     // dEnergy/du_z = sum_j (B_j(u_y) * dEnergy/du_z(y_j, u_z_interp))
        dot(py, (float4)(fe0.x, fe1.x, fe2.x, fe3.x))      // Energy       = sum_j (B_j(u_y) * Energy(y_j, u_z_interp))
    );
}

inline float4 fe3d_pbc_comb_tex(
    const float3 u,            // normalized coordinates (fractional voxel coordinates)
    const int3   n,            // dimensions of the B-spline grid (texture dimensions)
    __read_only image3d_t img, // 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
    const float4 PLQH          // coefficients to combine the 4 potential components
) {
    float4 bx,by,bz, dx,dy,dz;
    {  int ix = (int)u.x; const float tx=u.x - ix;  bx=basis(tx);  dx=dbasis(tx); }
    {  int iy = (int)u.y; const float ty=u.y - iy;  by=basis(ty);  dy=dbasis(ty); }
    {  int iz = (int)u.z; const float tz=u.z - iz;  bz=basis(tz);  dz=dbasis(tz); }
    const float3 dc = (float3){ 1.0f/n.x, 1.0f/n.y, 1.0f/n.z };
    const float3 uc = u*dc;
    // Interpolate along YZ for 4 different X planes
    // E#.x = dE/duy, E#.y = dE/duz, E#.z = E, all at (qx.#, u_y_interp, u_z_interp)
    const float3 E1 = fe2d_comb_tex(img, (float3){uc.x           , uc.y, uc.z}, PLQH, bz, dz, by, dy, dc);
    const float3 E2 = fe2d_comb_tex(img, (float3){uc.x+dc.x      , uc.y, uc.z}, PLQH, bz, dz, by, dy, dc);
    const float3 E3 = fe2d_comb_tex(img, (float3){uc.x+dc.x*2.0f , uc.y, uc.z}, PLQH, bz, dz, by, dy, dc);
    const float3 E4 = fe2d_comb_tex(img, (float3){uc.x+dc.x*3.0f , uc.y, uc.z}, PLQH, bz, dz, by, dy, dc);
    // Interpolate along X using results from 2D YZ-interpolation
    // Result is (dE/dux, dE/duy, dE/duz, E_total)
    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),     // dEnergy/du_x = sum_i (B'_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),     // dEnergy/du_y = sum_i (B_i(u_x) * dEnergy/du_y(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),     // dEnergy/du_z = sum_i (B_i(u_x) * dEnergy/du_z(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))      // Energy       = sum_i (B_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
    );
}


// ======================================================================
//                           getNonBond_GridFF_Bspline_tex()
// ======================================================================
// Calculate non-bonded forces on atoms (including both node atoms and capping atoms), considering periodic boundary conditions
// It calculates Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system,
// and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by tricubic B-spline texture interpolation
// It can be run in parallel for multiple systems.
// NOTE: This version uses a 3D texture for the B-spline GridFF.
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF_Bspline_tex( // Renamed kernel to distinguish from buffer version
    const int4 nDOFs,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  apos,         // 2 // positions of atoms
    __global float4*  aforce,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW_cutoff_factor_for_LJ, alphaMorse, Q_atom, H_bond_params_unused)
    // GridFF - Using Texture
    __read_only image3d_t BsplinePLQH_tex, // 10 // Grid-Force-Field (GFF) data (Pauli,London,Coulomb,HBond) in a 3D texture (Renamed texture)
    const int4     grid_ns,         // 11 // dimensions of the grid (matches buffer code name)
    const float4   grid_invStep,    // 12 // inverse of grid cell dimensions
    const float4   grid_p0          // 13 // origin of the grid
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms 
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=nDOFs.x;         // number of atoms in the system
    const int nnode =nDOFs.y;         // number of nodes in the system
    const int nvec  =natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system
    
    const float4 REQKi = REQs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = apos    [iav].xyz;       // position of the atom
    float4 fe          = float4Zero;              // forces on the atom

    const int iS_DBG = 0;
    const int iG_DBG = 0;

    // =================== Non-Bonded interaction ( molecule-molecule )

    { // insulate nbff

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "GPU::getNonBond_GridFF_Bspline() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){ 
    //     printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); 
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );        
    //     }
    // }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom

    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = apos [i+i0v];      // load positions  of atoms into local memory
        LCLJS [iL] = REQs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory 
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )   
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK.x  +=REQKi.x;              // mixing of RvdW radii
                REQK.yz *=REQKi.yz;             // mixing of EvdW and Q
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom            
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){      
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ? 
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    } // insulate nbff

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory 

    // ========== Molecule-Grid interaction with GridFF using tricubic Bspline (Texture based) ==================

    //__local int4 xqs[4]; // Local memory for PBC index patterns for X
    //__local int4 yqs[4]; // Local memory for PBC index patterns for Y

    { // insulate gridff
        // Initialize local memory for PBC index patterns. Only first 8 work-items do this.
        // ToDo: We can perhaps remove this when using hardware periodic boundary conditions ( sampler with CLK_ADDRESS_REPEAT )
        //if      (iL<4){ xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
        //else if (iL<8){ yqs[iL-4]=make_inds_pbc(grid_ns.y,iL-4); };
        //barrier(CLK_LOCAL_MEM_FENCE); // Ensure local memory is populated

        // Coefficients for combining Pauli, London, Coulomb, H-bond from the grid field
        // Matches original calculation using GFFParams.y (alphaMorse) and REQKi
        const float alphaMorse = GFFParams.y;
        const float ej = exp( alphaMorse * REQKi.x ); // REQKi.x is RvdW of atom_i
        const float4 PLQH = (float4){
            ej*ej*REQKi.y,                   // Pauli coeff:   EvdW_i * exp(2 * alphaMorse * RvdW_i)
            ej*   REQKi.y,                   // London coeff:  EvdW_i * exp(    alphaMorse * RvdW_i)
            REQKi.z,                         // Coulomb coeff: Q_i
            0.0f                             // H-bond coeff (assuming zeroed out)
        };
        //const float4 PLQH = (float4){0.0f,0.0f,0.0f,1.0f};

        // Calculate normalized coordinates 'u' for B-spline interpolation
        const float3 u  = (posi - grid_p0.xyz) * grid_invStep.xyz;
        //const float3 uc = (posi - grid_p0.xyz) * grid_invStep.xyz * (float3){ 1.f/(float)grid_ns.x, 1.f/(float)grid_ns.y, 1.f/(float)grid_ns.z};
        // Perform 3D B-spline interpolation using texture
        // fg contains (dE/dux, dE/duy, dE/duz, Energy)
        float4 fg = fe3d_pbc_comb_tex(u, grid_ns.xyz, BsplinePLQH_tex, PLQH);

        //float4 fg = fe3d_pbc_comb_tex(u, uc, BsplinePLQH_tex, PLQH);
        
        //fg.xyz *= -grid_invStep.xyz; // dux/dx = grid_invStep.x, etc.   Fx = -dE/dx = - (dE/dux) * (dux/dx)

        //float4 fg = read_imagef(BsplinePLQH_tex, sampler_bspline, (float4)( u.x, u.y, u.z, 0.0));
        //float4 fg = read_imagef(BsplinePLQH_tex, sampler_bspline, (int4)( (int)(u.x+0.5f), (int)(u.y+0.5f), (int)(u.z+0.5f), 0));

        //float4 fg = read_imagef(BsplinePLQH_tex, sampler_bspline, (float4)( posi.x, posi.y, posi.z, 0.0) );
        //float4 fg = read_imagef(BsplinePLQH_tex, sampler_bspline, (float4)( posi.x+0.5f, posi.y+0.5f, posi.z+0.5f, 0.0) );

        //float4 fg = read_imagef(BsplinePLQH_tex, sampler_norm_rep_near, (float4)( uc.x+0.000001, uc.y+0.000001, uc.z+0.000001, 0.f ));

        fe += fg;

        // fes[iG] = fe; // If you have a separate energy buffer

        //if((iG==iG_DBG)&&(iS==iS_DBG)){   printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8f,%16.8f,%16.8f,%16.8f) grid_p0(%16.8f,%16.8f,%16.8f) grid_invStep(%16.8f,%16.8f,%16.8f) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w, grid_p0.x,grid_p0.y,grid_p0.z, grid_invStep.x,grid_invStep.y,grid_invStep.z );   }
        //if((iG==iG_DBG)&&(iS==iS_DBG)){   printf( "GPU[%i] PLQH(%16.8e,%16.8e,%16.8e,%16.8e) apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8f,%16.8f,%16.8f,%16.8f) grid_ns(%3i,%3i,%3i) grid_invStep(%16.8f,%16.8f,%16.8f) \n", iG, PLQH.x,PLQH.y,PLQH.z,PLQH.w, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w, grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x,grid_invStep.y,grid_invStep.z );   }
        //if((iG==iG_DBG)&&(iS==iS_DBG)){   printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8f,%16.8f,%16.8f,%16.8f) grid_ns(%3i,%3i,%3i) grid_invStep(%16.8f,%16.8f,%16.8f) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w, grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x,grid_invStep.y,grid_invStep.z );   }
        
        //if((iG==iG_DBG)&&(iS==iS_DBG)){   printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8f,%16.8f,%16.8f,%16.8f) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w );   }
        // if((iG==iG_DBG)&&(iS==iS_DBG)){   
        //     //printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) fg(%16.8e,%16.8e,%16.8e,%16.8e) \n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, fg.x,fg.y,fg.z,fg.w );   
        //     //float4 feg = read_imagef(BsplinePLQH_tex, sampler_bspline, (int4)(u.x, -5,     u.z, 0));
        //     printf("GPU[%i] apos(%16.8f,%16.8f,%16.8f) u(%16.8f,%16.8f,%16.8f) feg(%16.8e,%16.8e,%16.8e,%16.8e)\n", iG, posi.x,posi.y,posi.z, u.x,u.y,u.z, feg.x,feg.y,feg.z,feg.w);
        // }
    }  // insulate gridff

    //if((iG==iG_DBG)&&(iS==iS_DBG)){   printf( "GPU[%i] apos(%16.8f,%16.8f,%16.8f) force(%16.8f,%16.8f,%16.8f,%16.8f) \n", iG, posi.x,posi.y,posi.z, fe.x,fe.y,fe.z,fe.w );}

    // Store the total force and energy for this atom
    aforce[iav] = fe;
    // Use forces[iav] += fe; if forces buffer accumulates from multiple kernels
}

__kernel void sampleGrid_tex(
    const int4 nDOFs,                      // 1
    __global float4*  apos,                // 2
    __global float4*  aforce,              // 3
    __global float4*  REQs,                // 4
    const float4  GFFParams,               // 5
    __read_only image3d_t BsplinePLQH_tex, // 6
    const float4   grid_invStep,           // 7
    const int4     grid_ns,                // 8
    const float4   grid_p0,                // 9
    const float4   MDparams,               // 10
    const int      nstep                   // 11
){
    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int np = nDOFs.x;

    float3 dz = (float3){ 0.0f, 0.0f, 0.1f };

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    // const float4 REQ        = REQs [iG];
    // const float  R2damp     = GFFParams.x*GFFParams.x;
    // const float  alphaMorse = GFFParams.y;

    float dt    =      MDparams.x;
    float cdamp = 1.0f-MDparams.y;

    // if(iG==0){ 
    //     printf( "GPU::sampleGrid_tex() grid_ns(%i,%i,%i) grid_invStep(%g,%g,%g) grid_p0(%g,%g,%g) \n", grid_ns.x,grid_ns.y,grid_ns.z, grid_invStep.x, grid_invStep.y, grid_invStep.z ); 
    //     printf( "GPU::sampleGrid_tex() np=%i nstep=%i dt=%g cdamp=%g grid_ns(%i,%i,%i,%i) \n", np, nstep, dt, cdamp, grid_ns.x,grid_ns.y,grid_ns.z,grid_ns.w ); 
    // }

    const float4 PLQH       = (float4){0.0f,0.0f,1.0f,0.0f};

    float3 posi  = apos[iG].xyz;
    float4 fe    = (float4){0.0f,0.0f,0.0f,0.0f};
    float3 vel   = (float3){0.0f,0.0f,0.0f};

    for(int i=0; i<nstep; i++){
        float3 u  = (posi - grid_p0.xyz) * grid_invStep.xyz;
        u        +=(float3){0.5f,0.5f,0.5f};
        fe        = fe3d_pbc_comb_tex(u, grid_ns.xyz, BsplinePLQH_tex, PLQH);
        vel       *= cdamp;
        vel       += fe.xyz*dt;
        posi      += vel   *dt;
    }

    aforce[iG]              = fe;

    //aforce[iG] = read_imagef( BsplinePLQH_tex, sampler_bspline, (int4)( (int)posi.x, (int)posi.y, (int)posi.z, 0 ) );
    //aforce[iG] = read_imagef( BsplinePLQH_tex, sampler_bspline, (float4)( posi.x, posi.y, posi.z, 0.0f ) );    

    // if(iG==0){ 
    //     printf( "GPU::sampleGridFF_tex() np=%i R2damp=%g aMorse=%g \n", np, R2damp, alphaMorse ); 
    //     for(int i=0; i<np; i++){
    //         float4 pi = apos[i];
    //         float4 fe = read_imagef( BsplinePLQH_tex, sampler_bspline, (float4)( pi.x, pi.y, pi.z, 0.0f ) );
    //         printf( "GPU::sampleGridFF_tex() %i p(%16.8f,%16.8f,%16.8f) force(%16.8f,%16.8f,%16.8f|%16.8f) \n", i, pi.x,pi.y,pi.z, fe.x,fe.y,fe.z,fe.w ); 
    //     }; 
    // }
// NOTE: https://registry.khronos.org/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
// CLK_ADDRESS_REPEAT - out-of-range image coordinates are wrapped to the valid range. This address mode can only be used with normalized coordinates. If normalized coordinates are not used, this addressing mode may generate image coordinates that are undefined.
}