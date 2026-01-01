// =================================================================================
// PHYSICS CONSTANTS & MACROS
// =================================================================================
#define TILE_W 8
#define TILE_H 8
#define HALO   1
#define LOCAL_W (TILE_W + 2 * HALO)
#define LOCAL_H (TILE_H + 2 * HALO)
#define LOCAL_SIZE (LOCAL_W * LOCAL_H)

// Material 1: Liquid Hydrogen
// Real density ~71 kg/m^3. Gamma ~1.4.
#define GAMMA_1 1.4f
#define PI_1    0.0f 

// Material 2: Uranium
// Real density ~19050 kg/m^3. Gamma ~3.0. P_inf ~40 GPa.
#define GAMMA_2 3.0f  
#define PI_2    40.0e9f 

#define EPS 1e-6f
#define MIN_RHO 1e-4f
#define MIN_E   1e-4f

typedef struct {
    float rho_a1; 
    float rho_a2; 
    float ru;     
    float rv;     
    float E;      
    float phi;    
} Conserved;

typedef struct {
    float p;      
    float u;      
    float v;      
    float c;      
} Primitive;

// =================================================================================
// HELPER FUNCTIONS
// =================================================================================

inline float heaviside(float phi, float dx) {
    float epsilon = 1.5f * dx;
    float x = (phi + epsilon) / (2.0f * epsilon);
    x = clamp(x, 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
}

inline float smoothed_sign(float phi, float dx) {
    return phi / sqrt(phi*phi + dx*dx);
}

// =================================================================================
// ROBUST EOS
// =================================================================================
inline Primitive get_primitive(Conserved U, float dx) {
    Primitive P;
    
    // 1. Recover Density and handle vacuum protection
    float rho = U.rho_a1 + U.rho_a2;
    if (rho < MIN_RHO) {
        P.u = 0.0f;
        P.v = 0.0f;
        P.p = MIN_E;
        P.c = 340.0f; 
        return P;
    }

    // 2. Recover Velocities
    P.u = U.ru / rho;
    P.v = U.rv / rho;
    float ke = 0.5f * (P.u*P.u + P.v*P.v);

    // 3. Recover Internal Energy with floor
    float internal_e = U.E - rho * ke;
    if (internal_e < MIN_E) internal_e = MIN_E;

    // 4. Mixture Equation of State (using Level-Set for Volume Fraction)
    float alpha2 = heaviside(U.phi, dx);
    float alpha1 = 1.0f - alpha2;

    float g1_m1 = GAMMA_1 - 1.0f;
    float g2_m1 = GAMMA_2 - 1.0f;

    float denom = (alpha1 / g1_m1) + (alpha2 / g2_m1) + EPS;
    float pinf_term = (alpha1 * GAMMA_1 * PI_1 / g1_m1) + (alpha2 * GAMMA_2 * PI_2 / g2_m1);

    P.p = (internal_e - pinf_term) / denom;
    if (P.p < 1.0f) P.p = 1.0f; // Minimal pressure

    // 5. Speed of Sound (Robust mixture sound speed)
    // Using the frozen sound speed definition for 5-equation model
    float gamma_avg_num = (alpha1 * GAMMA_1 / g1_m1) + (alpha2 * GAMMA_2 / g2_m1);
    float c_sq = (gamma_avg_num / denom) * (P.p + pinf_term / gamma_avg_num) / (rho + EPS);
    
    P.c = sqrt(max(c_sq, 1.0f));

    return P;
}

// =================================================================================
// FLUX & UPDATE
// =================================================================================

inline Conserved compute_flux(Conserved UL, Conserved UR, float dx, float dt, float nx, float ny) {
    Primitive PL = get_primitive(UL, dx);
    Primitive PR = get_primitive(UR, dx);

    float un_L = PL.u * nx + PL.v * ny;
    float un_R = PR.u * nx + PR.v * ny;

    float S_L = fabs(un_L) + PL.c;
    float S_R = fabs(un_R) + PR.c;
    
    // SAFETY: Add 10% to wave speed to ensure stability (Entropy fix)
    float S_max = max(S_L, S_R) * 1.1f;

    Conserved FL, FR;

    // Flux L
    FL.rho_a1 = UL.rho_a1 * un_L;
    FL.rho_a2 = UL.rho_a2 * un_L;
    FL.ru     = UL.ru * un_L + PL.p * nx;
    FL.rv     = UL.rv * un_L + PL.p * ny;
    FL.E      = (UL.E + PL.p) * un_L;
    FL.phi    = UL.phi * un_L;

    // Flux R
    FR.rho_a1 = UR.rho_a1 * un_R;
    FR.rho_a2 = UR.rho_a2 * un_R;
    FR.ru     = UR.ru * un_R + PR.p * nx;
    FR.rv     = UR.rv * un_R + PR.p * ny;
    FR.E      = (UR.E + PR.p) * un_R;
    FR.phi    = UR.phi * un_R;

    Conserved F;
    F.rho_a1 = 0.5f * (FL.rho_a1 + FR.rho_a1 - S_max * (UR.rho_a1 - UL.rho_a1));
    F.rho_a2 = 0.5f * (FL.rho_a2 + FR.rho_a2 - S_max * (UR.rho_a2 - UL.rho_a2));
    F.ru     = 0.5f * (FL.ru     + FR.ru     - S_max * (UR.ru     - UL.ru));
    F.rv     = 0.5f * (FL.rv     + FR.rv     - S_max * (UR.rv     - UL.rv));
    F.E      = 0.5f * (FL.E      + FR.E      - S_max * (UR.E      - UL.E));
    F.phi    = 0.5f * (FL.phi    + FR.phi    - S_max * (UR.phi    - UL.phi));

    return F;
}

__kernel void update_fluid(
    __global const float4* restrict in_U_vec, 
    __global const float2* restrict in_E_vec,
    __global float4* restrict out_U_vec,
    __global float2* restrict out_E_vec,
    __global float*  restrict out_p,  // [Diagnostic] Pressure    , no need for calculation
    __global float*  restrict out_c,  // [Diagnostic] Sound speed , no need for calculation
    int width,
    int height,
    float dx,
    float dt
) {
    int lx = get_local_id(0); 
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;

    __local float l_rho_a1[LOCAL_SIZE];
    __local float l_rho_a2[LOCAL_SIZE];
    __local float l_ru[LOCAL_SIZE];
    __local float l_rv[LOCAL_SIZE];
    __local float l_E[LOCAL_SIZE];
    __local float l_phi[LOCAL_SIZE];

    int total_items = LOCAL_W * LOCAL_H; 
    int threads = get_local_size(0); 

    // LOAD
    for (int t = lx; t < total_items; t += threads) {
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;
        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;
        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        int g_idx = gy * width + gx;

        float4 vecA = in_U_vec[g_idx];
        float2 vecB = in_E_vec[g_idx];
        l_rho_a1[t] = vecA.x;
        l_rho_a2[t] = vecA.y;
        l_ru[t]     = vecA.z;
        l_rv[t]     = vecA.w;
        l_E[t]      = vecB.x;
        l_phi[t]    = vecB.y;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // COMPUTE
    int pixels_to_process = TILE_W * TILE_H; 
    
    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;
        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;

        Conserved Uc;
        Uc.rho_a1 = l_rho_a1[c_idx];
        Uc.rho_a2 = l_rho_a2[c_idx];
        Uc.ru     = l_ru[c_idx];
        Uc.rv     = l_rv[c_idx];
        Uc.E      = l_E[c_idx];
        Uc.phi    = l_phi[c_idx];

        // Vacuum momentum kill
        if (Uc.rho_a1 + Uc.rho_a2 < MIN_RHO) {
             Uc.ru = 0.0f; Uc.rv = 0.0f;
        }

        #define GET_U(idx) (Conserved){l_rho_a1[idx], l_rho_a2[idx], l_ru[idx], l_rv[idx], l_E[idx], l_phi[idx]}
        Conserved U_left  = GET_U(cy * LOCAL_W + (cx - 1));
        Conserved U_right = GET_U(cy * LOCAL_W + (cx + 1));
        Conserved U_down  = GET_U((cy - 1) * LOCAL_W + cx);
        Conserved U_up    = GET_U((cy + 1) * LOCAL_W + cx);

        Conserved F_xm = compute_flux(U_left, Uc, dx, dt, 1.0f, 0.0f);  
        Conserved F_xp = compute_flux(Uc, U_right, dx, dt, 1.0f, 0.0f); 
        Conserved F_ym = compute_flux(U_down, Uc, dx, dt, 0.0f, 1.0f);  
        Conserved F_yp = compute_flux(Uc, U_up, dx, dt, 0.0f, 1.0f);    

        float dtdx = dt / dx; 

        Conserved U_new;
        U_new.rho_a1 = Uc.rho_a1 - dtdx * (F_xp.rho_a1 - F_xm.rho_a1 + F_yp.rho_a1 - F_ym.rho_a1);
        U_new.rho_a2 = Uc.rho_a2 - dtdx * (F_xp.rho_a2 - F_xm.rho_a2 + F_yp.rho_a2 - F_ym.rho_a2);
        U_new.ru     = Uc.ru     - dtdx * (F_xp.ru     - F_xm.ru     + F_yp.ru     - F_ym.ru);
        U_new.rv     = Uc.rv     - dtdx * (F_xp.rv     - F_xm.rv     + F_yp.rv     - F_ym.rv);
        U_new.E      = Uc.E      - dtdx * (F_xp.E      - F_xm.E      + F_yp.E      - F_ym.E);
        // 4. Update Level Set (phi) - NON-CONSERVATIVE ADVECTION
        // phi is a signed distance field, not a density. Conservative update is wrong.
        // dphi/dt + u*dphi/dx + v*dphi/dy = 0
        Primitive Pc = get_primitive(Uc, dx);
        float grad_phi_x = (U_right.phi - U_left.phi) / (2.0f * dx);
        float grad_phi_y = (U_up.phi    - U_down.phi) / (2.0f * dx);
        U_new.phi = Uc.phi - dt * (Pc.u * grad_phi_x + Pc.v * grad_phi_y);
        
        // --- POSITIVITY FLOORS (Prevent NaNs) ---
        // 1. Partial Density Floors
        if (U_new.rho_a1 < 0.0f) U_new.rho_a1 = MIN_RHO * 0.5f;
        if (U_new.rho_a2 < 0.0f) U_new.rho_a2 = MIN_RHO * 0.5f;
        float rho_tot = U_new.rho_a1 + U_new.rho_a2;
        if (rho_tot < MIN_RHO) {
            U_new.rho_a1 = MIN_RHO * 0.5f;
            U_new.rho_a2 = MIN_RHO * 0.5f;
            rho_tot = MIN_RHO;
        }

        // 2. Kinetic Energy / Momentum clamping
        float vel_sq = (U_new.ru * U_new.ru + U_new.rv * U_new.rv) / (rho_tot * rho_tot + EPS);
        float ke_new = 0.5f * rho_tot * vel_sq;

        // 3. Energy Floor (ensure internal energy is positive)
        // We need internal_e = E - ke > MIN_E
        if (U_new.E < ke_new + MIN_E) {
            U_new.E = ke_new + MIN_E;
        }

        // 4. Vacuum Momentum Kill (Zero velocities in low-density cells)
        if (rho_tot < MIN_RHO * 1.1f) {
             U_new.ru = 0.0f; 
             U_new.rv = 0.0f;
        }

        // Output with bounds checks
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;

        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            out_U_vec[out_idx] = (float4)(U_new.rho_a1, U_new.rho_a2, U_new.ru, U_new.rv);
            out_E_vec[out_idx] = (float2)(U_new.E, U_new.phi);
            Primitive P_new = get_primitive(U_new, dx);
            out_p[out_idx] = P_new.p;
            out_c[out_idx] = P_new.c;
        }
    }
}

__kernel void redistance_phi(
    __global const float2* restrict in_E_vec,  
    __global const float2* restrict orig_E_vec,
    __global float2* restrict out_E_vec,       
    int width,
    int height,
    float dx,
    float dtau 
) {
    int lx = get_local_id(0);
    int grp_x = get_group_id(0);
    int grp_y = get_group_id(1);
    int tile_origin_x = grp_x * TILE_W; 
    int tile_origin_y = grp_y * TILE_H;
    __local float l_phi[LOCAL_SIZE];
    int total_items = LOCAL_W * LOCAL_H;
    int threads = get_local_size(0);
    for (int t = lx; t < total_items; t += threads) {
        int ly_loc = t / LOCAL_W;
        int lx_loc = t % LOCAL_W;
        int gx = tile_origin_x + lx_loc - HALO;
        int gy = tile_origin_y + ly_loc - HALO;
        gx = max(0, min(gx, width - 1));
        gy = max(0, min(gy, height - 1));
        l_phi[t] = in_E_vec[gy * width + gx].y;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    int pixels_to_process = TILE_W * TILE_H;
    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;
        int cx = tx + HALO;
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;
        float phi_c = l_phi[c_idx];
        float phi_W = l_phi[cy * LOCAL_W + (cx - 1)]; 
        float phi_E = l_phi[cy * LOCAL_W + (cx + 1)]; 
        float phi_S = l_phi[(cy - 1) * LOCAL_W + cx]; 
        float phi_N = l_phi[(cy + 1) * LOCAL_W + cx]; 
        float dx_m = (phi_c - phi_W) / dx;
        float dx_p = (phi_E - phi_c) / dx;
        float dy_m = (phi_c - phi_S) / dx;
        float dy_p = (phi_N - phi_c) / dx;
        int gx_out = tile_origin_x + tx;
        int gy_out = tile_origin_y + ty;
        float phi_0 = orig_E_vec[gy_out * width + gx_out].y;
        float S = smoothed_sign(phi_0, dx);
        float grad_phi_sq = 0.0f;
        if (phi_0 > 0.0f) {
            float dx_term = max(max(0.0f, -dx_p), max(0.0f, dx_m)); 
            float dy_term = max(max(0.0f, -dy_p), max(0.0f, dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        } else {
            float dx_term = max(max(0.0f, dx_p), max(0.0f, -dx_m));
            float dy_term = max(max(0.0f, dy_p), max(0.0f, -dy_m));
            grad_phi_sq = dx_term*dx_term + dy_term*dy_term;
        }
        float grad_mag = sqrt(grad_phi_sq);
        float phi_new = phi_c - dtau * S * (grad_mag - 1.0f);
        if(gx_out < width && gy_out < height) {
            int out_idx = gy_out * width + gx_out;
            float E_val = in_E_vec[out_idx].x; 
            out_E_vec[out_idx] = (float2)(E_val, phi_new);
        }
    }
}
