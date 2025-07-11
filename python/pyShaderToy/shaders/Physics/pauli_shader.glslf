// Shadertoy GLSL shader for ppafm simplified current calculation

// --- Helper Functions ---


// #define NSINGLE   3
// #define NSTATE    8

//#define NSINGLE   4
//#define NSTATE    16

#define NSINGLE   6
#define NSTATE    64

// Basic vector operations
vec3  vec3_sub(vec3 a, vec3 b) { return a - b; }
float vec3_norm(vec3 v) { return length(v); }
float vec3_norm2(vec3 v) { return dot(v,v); }

// Bit manipulation
bool site_in_state(int site, int state) { return ((state >> site) & 1) != 0; }

int count_electrons(int state) {
    int count = 0;
    for (int i = 0; i < NSINGLE; ++i) {
        if (site_in_state(i, state)) count++;
    }
    return count;
}

// Constants
const float COULOMB_CONST = 14.3996448915; // [eV A]




/*
// Emultipole function
float Emultipole(const vec3 d, int order, const float cs[10]) {
    float ir2 = 1.0/vec3_norm2(d);
    float E = cs[0];
    if(order>0) E += ir2*(cs[1]*d.x + cs[2]*d.y + cs[3]*d.z);
    if(order>1) E += ir2*ir2*((cs[4]*d.x + cs[9]*d.y)*d.x +
                             (cs[5]*d.y + cs[7]*d.z)*d.y +
                             (cs[6]*d.z + cs[8]*d.x)*d.z);
    return sqrt(ir2)*E;
}

// evalMultipoleMirror function
float evalMultipoleMirror(vec3 pTip, const vec3 pSite, float VBias, float Rtip, vec2 zV, int order, const float cs[10], float E0, mat3 rotSite, bool bMirror, bool bRamp) {
    float zV0 = zV.x;
    float zVd = zV.y;
    float orig_z = pTip.z;
    float zV1 = orig_z + zVd;

    vec3 pTipMirror = pTip;
    pTipMirror.z = 2.0*zV0 - orig_z;

    pTip = pTip - pSite;
    pTipMirror = pTipMirror - pSite;
    
    pTip = rotSite * pTip;
    pTipMirror = rotSite * pTipMirror;

    float E = Emultipole(pTip, order, cs);
    if(bMirror) E -= Emultipole(pTipMirror, order, cs);
    float VR = VBias * Rtip;
    E *= VR;

    if(bRamp) {
        float ramp = (pSite.z - zV0) / (zV1 - zV0);
        ramp = clamp(ramp, 0.0, 1.0);
        if(pSite.z < zV0) ramp = 0.0;
        float V_lin = VBias * ramp;
        float E_lin = cs[0] * V_lin;
        E += E_lin;
    }
    
    return E + E0;
}
*/

float evalEtip_simple(vec3 pTip, const vec3 pSite, float VBias, float Rtip ) {
    float r = vec3_norm(pTip-pSite);
    float E = VBias * Rtip / r;
    return E;
}

// evalSitesTipsTunneling function
float evalSitesTipsTunneling(const vec3 pTip, const vec3 pSite, float beta, float Amp) {
    vec3 d = pTip - pSite;
    float r = length(d);
    return Amp * exp(-beta * r);
}


// PauliSolver functions

float calculate_state_energy_simple(int state, const float Esingle[NSINGLE], float W) {
    float E = 0.0;
    for(int i=0; i<NSINGLE; i++) {
        if(site_in_state(i, state)){ 
            E += Esingle[i];
            for(int j=0; j<NSINGLE; j++) {
                if(j==i) continue;
                if(site_in_state(j, state)) E += W;
            }
        }
    }
    return E;
}


float calculate_state_energy_dist(int state, const float Esingle[NSINGLE], const vec3 pSites[NSINGLE], float W) {
    float E = 0.0;
    for(int i=0; i<NSINGLE; i++) {
        if(site_in_state(i, state)){ 
            E += Esingle[i];
            for(int j=0; j<NSINGLE; j++) {
                if(j==i) continue;
                if(site_in_state(j, state)) E += W*vec3_norm(pSites[i]-pSites[j]);
            }
        }
    }
    return E;
}



float solve_equilibrium_ground_state( inout float energies[NSTATE], inout float probabilities[NSTATE], inout float Qs[NSINGLE]) {
    int idx_min = 0;
    float e_min = energies[0];
    for(int i=1; i<NSTATE; i++) {
        probabilities[i] = 0.0;
        if(energies[i] < e_min) { e_min = energies[i]; idx_min = i; }

    }
    float Qtot = 0.0;
    for(int i=0; i<NSINGLE; i++) {
        Qs[i] = site_in_state(i, idx_min) ? 1.0 : 0.0;
        Qtot += Qs[i];
    }
    probabilities[idx_min] = 1.0;
    return Qtot;
}

/*
void solve_equilibrium_boltzmann( float temperature, float chemical_potential, int nSingle, inout float energies[NSTATE], inout float probabilities[NSTATE]) {
    if(temperature <= 0.0) {
        solve_equilibrium_ground_state(energies, probabilities);
        return;
    }

    float exponents[8];
    float max_exp = -1e30;
    for(int i = 0; i < nstates; ++i) {
        int N = count_electrons(i, nSingle);
        float exponent = -(energies[i] - float(N)*chemical_potential)/temperature;
        exponents[i] = exponent;
        if(exponent > max_exp) max_exp = exponent;
    }

    float Z = 0.0;
    for(int i = 0; i < nstates; ++i) {
        float w = exp(exponents[i] - max_exp);
        Z += w;
    }

    if(!(Z > 0.0) || isinf(Z)) {
        solve_equilibrium_ground_state(energies, probabilities);
        return;
    }

    for(int i = 0; i < NSTATE; ++i) {
        probabilities[i] = exponents[i] / Z;
    }
}

float current_simple(int nSingle, const float probabilities[NSTATE], const float TLeads[3], bool bEmpty) {
    float occ[3] = float[3](0.0, 0.0, 0.0);
    
    for(int state = 0; state < NSTATE; ++state) {
        float p = probabilities[state];
        int mask = state;
        for(int s = 0; s < nSingle; ++s) {
            if((mask & 1) != 0) occ[s] += p;
            mask >>= 1;
        }
    }
    
    float current = 0.0;
    for(int s = 0; s < nSingle; ++s) {
        float p = occ[s];
        if(bEmpty) p = 1.0-p;
        current += TLeads[s] * p;
    }
    return current;
}
*/

// --- Hardcoded Parameters (Example values - will be replaced with actual data) ---
// const int nSingle = 3; // Number of single-particle states (e.g., sites)
// const int nstates = 8; // Total many-body states (2^nSingle)
const int order = 1; // Multipole order

#define Rmol  5.0
#define Rmol2 7.0

// Example site positions (replace with actual data)
// const vec3 pSites[NSINGLE] = vec3[NSINGLE](
//     vec3(-0.5*Rmol, -0.866*Rmol , 0.0),
//     vec3( 1.0*Rmol,      0.0    , 0.0),
//     vec3(-0.5*Rmol,  0.866*Rmol , 0.0)
// );
// float Es0[NSINGLE] = float[NSINGLE](-0.2, -0.2, -0.2);

// const vec3 pSites[NSINGLE] = vec3[NSINGLE](
//     vec3(-1.0*Rmol, -1.0*Rmol , 0.0),
//     vec3( 1.0*Rmol, -1.0*Rmol , 0.0),
//     vec3(-1.0*Rmol,  1.0*Rmol , 0.0),
//     vec3( 1.0*Rmol,  1.0*Rmol , 0.0)
// );
// float Es0[NSINGLE] = float[NSINGLE](-0.2, -0.2, -0.2, -0.2);

const vec3 pSites[NSINGLE] = vec3[NSINGLE](
    vec3(-0.5*Rmol,  -0.866*Rmol  , 0.0),
    vec3( 1.0*Rmol,      0.0      , 0.0),
    vec3(-0.5*Rmol,   0.866*Rmol  , 0.0),
    vec3(+0.5*Rmol2,  0.866*Rmol2 , 0.0),
    vec3(-1.0*Rmol2,      0.0     , 0.0),
    vec3(+0.5*Rmol2, -0.866*Rmol2 , 0.0)
);
float Es0[NSINGLE] = float[NSINGLE]( -0.2, -0.2, -0.2,   -0.2,-0.2, -0.2 );


// Example site rotations (identity matrices for now)
// const mat3 rotSites[3] = mat3[3](
//     mat3(1.0),
//     mat3(1.0),
//     mat3(1.0)
// );

// Example multipole coefficients (cs[0] is monopole, cs[1-3] dipole, cs[4-9] quadrupole)
// const float cs[10] = float[10](
//     1.0, 0.0, 0.0, 0.0, // Monopole, Dipole
//     0.0, 0.0, 0.0, 0.0, 0.0, 0.0  // Quadrupole
// );



// Example parameters
const float Rtip  = 2.0;
const float VBias = 0.7;
const float E0    = 0.0;
const float beta  = 1.0;
const float Amp   = 1.0;
const float W     = 0.02;
const float mu    = 0.0;
const float T     = 0.025;
const vec2 zV     = vec2(0.0, 0.0);
const int which_solver = -1; 


void mainImage( out vec4 fragColor, in vec2 fragCoord ){

    // --- Tip Position from fragCoord ---
    vec3 pTip = vec3(fragCoord.x / iResolution.x - 0.5, fragCoord.y / iResolution.y -0.5, 0.0);
    pTip.xy*=40.0;

    

    float Es[NSINGLE];
    float Ts[NSINGLE];
    float Qs[NSINGLE];
    float energies     [NSTATE];
    float probabilities[NSTATE];
 
    // 1. Calculate Es and Ts for each site
    float Emax = -1e+10;
    float Tmax = -1e+10;
    for (int i = 0; i < NSINGLE; ++i) {
        //float Ei = evalMultipoleMirror(pTip, pSites[i], VBias, Rtip, zV, order, cs, E0, rotSites[i], true, true);
        float Ei = Es0[i]+evalEtip_simple(pTip, pSites[i], VBias, Rtip);
        float Ti = evalSitesTipsTunneling(pTip, pSites[i], beta, Amp );
        //float Ti = 0.1/vec3_norm(pTip-pSites[i]);
        Es[i] = Ei; 
        Ts[i] = Ti;
        Emax = max(Emax, Ei);
        Tmax = max(Tmax, Ti);
    }

    float EmaxManyBody = -1e+10;
    float EminManyBody = 1e+10;
    for (int state=0; state<NSTATE; state++) {
        //float E = calculate_state_energy_simple(state, Es, W);
        float E = calculate_state_energy_dist(state, Es, pSites, W*0.1);
        energies[state] = E;
        EmaxManyBody = max(EmaxManyBody, E);
        EminManyBody = min(EminManyBody, E);
    }

    float Qsum = solve_equilibrium_ground_state( energies, probabilities, Qs);

    // int idx_min = 0;
    // float e_min = energies[0];
    // for(int i=1; i<nstates; i++) {
    //     probabilities[i] = 0.0;
    //     if(energies[i] < e_min) { e_min = energies[i]; idx_min = i; }

    // }
    // for(int i=0; i<NSINGLE; i++) {
    //     Qs[i] = site_in_state(i, idx_min) ? 1.0 : 0.0;
    // }
    // probabilities[idx_min] = 1.0;


    // 3. Solve for Probabilities
    // if (which_solver == -1) {
    //     solve_equilibrium_ground_state(nstates, energies, probabilities);
    // } else if (which_solver == -2) {
    //     solve_equilibrium_boltzmann(nstates, T, mu, NSINGLE, energies, probabilities);
    // } else {
    //     for(int i=0; i<nstates; ++i) {
    //         probabilities[i] = 0.0;
    //     }
    // }

    /*
    // 4. Calculate Current
    float current = current_simple(NSINGLE, nstates, probabilities, Ts, true);

    // --- Output Visualization ---
    float display_current = current * 1000.0;
    float normalized_current = (display_current + 1.0) / 2.0;
    normalized_current = clamp(normalized_current, 0.0, 1.0);
    */

    //vec3 color = mix(vec3(0.0, 0.0, 1.0), vec3(1.0, 0.0, 0.0), normalized_current);
    //vec3 color = vec3(Tmax, 0.0, 0.0);
    //vec3 color = vec3(Ts[0], Ts[1], Ts[2]);
    //vec3 color = vec3(Es[0], Es[1], Es[2])*1.0;
    //vec3 color = vec3(Emax, 0.0, 0.0);
    //vec3 color = vec3(s[0], Es[1], Es[2])*1.0;
    //vec3 color = vec3(Qs[0], Qs[1], Qs[2])*1.0;
    vec3 color = vec3( EmaxManyBody*0.0, (EminManyBody+0.5)*0.0, Qsum*0.3);
    //vec3 color = vec3( energies[0], energies[1], energies[2]);
    //vec3 color = pTip-pSites[1];

    //vec3 pTip = vec3(fragCoord.x / iResolution.x, fragCoord.y / iResolution.y, 5.0);
    fragColor = vec4(color, 1.0);
}
