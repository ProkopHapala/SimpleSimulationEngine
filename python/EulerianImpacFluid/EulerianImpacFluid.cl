/// === AUTO-DOC BEGIN ===
/// @file EulerianImpacFluid.cl
/// @brief OpenCL kernels for 2D Eulerian multi-material compressible flow with level-set interface tracking.
///
/// Implements a 5-equation diffuse-interface model for two-material impact simulations
/// (e.g. uranium projectile hitting liquid hydrogen). Uses a **stiffened-gas EOS** per
/// material with a level-set phi field to determine volume fractions — this is critical
/// because mass-fraction-derived alphas give wrong pressure at material interfaces.
///
/// **Key design decisions**:
/// - **Lax-Friedrichs (local Lax-Friedrichs/Rusanov) flux** — simple, robust for strong
///   shocks; avoids Riemann solver complexity at the cost of some diffusion.
/// - **Non-conservative advection of phi** — phi is a signed-distance field, not a
///   density. Conservative fluxing smears it and breaks volume-fraction consistency.
/// - **8x8 tiles with 2-cell halo in local memory** — cache-friendly 2D stencil;
///   64 threads per group cooperatively load the 12x12 halo region.
///   2-cell halo needed for MUSCL 2nd-order reconstruction (slopes use i-1..i+2).
/// - **Positivity floors** on partial densities and internal energy — prevent NaNs
///   in vacuum cells and at strong-shock interfaces where float32 precision is
///   insufficient for the stiffened-gas EOS with large PI_2 (40 GPa).
/// - **Redistancing kernel** reinitializes phi to a signed distance function using
///   Godunov-type Hamilton-Jacobi solver — prevents interface thickening over time.
///
/// **Buffer layout**: U packed as float4(rho_a1, rho_a2, ru, rv), E packed as float2(E, phi).
/// Ping-pong between two buffer pairs each step.
/// === AUTO-DOC END ===

// =================================================================================
// PHYSICS CONSTANTS & MACROS
// =================================================================================
#define TILE_W 8
#define TILE_H 8
#define HALO   2  ///< MUSCL needs 2-cell halo for slope computation (i-1, i, i+1, i+2)
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

/// @brief Smoothed (cubic Hermite) Heaviside — maps level-set phi to volume fraction alpha in [0,1].
/// Smearing width is 1.5*dx, matching the interface thickness used in EOS and initialization.
inline float heaviside(float phi, float dx) {
    float epsilon = 1.5f * dx;
    float x = (phi + epsilon) / (2.0f * epsilon);
    x = clamp(x, 0.0f, 1.0f);
    return x * x * (3.0f - 2.0f * x);
}

/// @brief Regularized sign function for redistancing — avoids discontinuity at phi=0.
inline float smoothed_sign(float phi, float dx) {
    return phi / sqrt(phi*phi + dx*dx);
}

// =================================================================================
// ROBUST EOS
// =================================================================================
/// @brief Convert conserved variables to primitives (p, u, v, c) via mixture stiffened-gas EOS.
/// Volume fractions derived from phi (not mass fractions) — essential for correct interface pressure.
/// Includes vacuum protection and internal-energy floor to handle low-density cells.
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

/// @brief Local Lax-Friedrichs (Rusanov) numerical flux between left/right states along normal (nx,ny).
/// Uses max wave speed with 10% entropy fix for stability at strong shocks.
inline Conserved compute_flux(Conserved UL, Conserved UR, float dx, float dt, float nx, float ny) {
    Primitive PL = get_primitive(UL, dx);
    Primitive PR = get_primitive(UR, dx);

    float un_L = PL.u * nx + PL.v * ny;
    float un_R = PR.u * nx + PR.v * ny;

    float S_L = fabs(un_L) + PL.c;
    float S_R = fabs(un_R) + PR.c;
    
    // SAFETY: Add 5% to wave speed to ensure stability (Entropy fix)
    float S_max = max(S_L, S_R) * 1.05f;

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

/// @brief Rusanov flux using pre-computed primitives — avoids redundant EOS calls.
inline Conserved compute_flux_from_primitive(Conserved UL, Conserved UR, Primitive PL, Primitive PR, float dx, float dt, float nx, float ny) {
    float un_L = PL.u * nx + PL.v * ny;
    float un_R = PR.u * nx + PR.v * ny;

    float S_L = fabs(un_L) + PL.c;
    float S_R = fabs(un_R) + PR.c;
    float S_max = max(S_L, S_R) * 1.05f;

    Conserved FL, FR;

    FL.rho_a1 = UL.rho_a1 * un_L;
    FL.rho_a2 = UL.rho_a2 * un_L;
    FL.ru     = UL.ru * un_L + PL.p * nx;
    FL.rv     = UL.rv * un_L + PL.p * ny;
    FL.E      = (UL.E + PL.p) * un_L;
    FL.phi    = UL.phi * un_L;

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

// =================================================================================
// MUSCL 2nd-ORDER RECONSTRUCTION
// =================================================================================
// References:
//   [1] van Leer, B. (1979). "Towards the Ultimate Conservative Difference Scheme V:
//       A Second-Order Sequel to Godunov's Method." J. Comput. Phys. 32, 101-136.
//   [2] Toro, E.F. (2009). "Riemann Solvers and Numerical Methods for Fluid Dynamics,"
//       3rd ed., Springer. Chapter 13-14 (MUSCL-Hancock, slope limiters).
//   [3] Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters for Hyperbolic
//       Conservation Laws." SIAM J. Numer. Anal. 21(5), 995-1011.
//
// MUSCL = Monotone Upstream-Centered Schemes for Conservation Laws.
// Instead of piecewise-constant cell values (Godunov 1st order), we reconstruct
// piecewise-linear profiles within each cell using limited slopes:
//
//   U(x) = U_i + sigma_i * (x - x_i)    for x in [x_{i-1/2}, x_{i+1/2}]
//
// where sigma_i is the limited slope. At face (i, i+1):
//   U_L = U_i     + 0.5 * sigma_i      (left state, extrapolated from cell i)
//   U_R = U_{i+1} - 0.5 * sigma_{i+1}  (right state, extrapolated from cell i+1)
//
// The minmod limiter ensures TVD (Total Variation Diminishing) property:
//   sigma_i = minmod( (U_i - U_{i-1}) , (U_{i+1} - U_i) )
//   minmod(a,b) = sign(a)*min(|a|,|b|) if a*b > 0, else 0
//
// This gives 2nd-order accuracy in smooth regions and degrades to 1st-order
// (sigma=0) near discontinuities — no spurious oscillations (Godunov's theorem).

/// @brief minmod slope limiter — TVD limiter that selects the smaller slope when signs agree, zero otherwise.
/// Ref [3] Sweby (1984). Returns 0 if a and b have opposite signs (discontinuity detected).
inline float minmod(float a, float b) {
    return (a * b > 0.0f) ? ((fabs(a) < fabs(b)) ? a : b) : 0.0f;
}

/// @brief MUSCL face reconstruction: compute limited slope and extrapolate to face.
/// Given cell-centered values U_{i-1}, U_i, U_{i+1}, returns the left face state at (i+1/2).
/// U_face = U_i + 0.5 * minmod(U_i - U_{i-1}, U_{i+1} - U_i)
/// Ref [1] van Leer (1979), Ref [2] Toro (2009) §13.2.
inline float muscl_left(float U_im1, float U_i, float U_ip1) {
    float slope = minmod(U_i - U_im1, U_ip1 - U_i);
    return U_i + 0.5f * slope;
}

/// @brief MUSCL face reconstruction: right state at face (i-1/2) from cell i.
/// U_face = U_i - 0.5 * minmod(U_i - U_{i-1}, U_{i+1} - U_i)
inline float muscl_right(float U_im1, float U_i, float U_ip1) {
    float slope = minmod(U_i - U_im1, U_ip1 - U_i);
    return U_i - 0.5f * slope;
}

/// @brief Reconstruct PRIMITIVE variables at a face using MUSCL with minmod limiter.
/// Reconstructs (rho, u, v, p, phi) instead of conserved variables — this is critical
/// for multi-material problems where conserved-variable reconstruction creates
/// inconsistent (rho, E) combinations at material interfaces, leading to huge
/// pressure overshoots. Primitive reconstruction ensures p is bounded by neighbors.
/// Ref [2] Toro (2009) §14.3 — primitive vs conserved reconstruction;
/// Ref [4] Saurel & Abgrall (1999). "A Simple Method for Compressible Multifluid Flows."
///      SIAM J. Sci. Comput. 21(3), 1115-1145. — primitive reconstruction for multi-material.
///
/// @param l_rho  Total density (l_rho_a1 + l_rho_a2) cached in local memory
/// @param l_p, l_u, l_v, l_phi  Cached primitive variables in local memory
/// @param idx_m2, idx_m1, idx_C, idx_p1  4-cell stencil indices in local memory
/// @param U_left, U_right  Output: reconstructed conserved states at the face
/// @param P_left, P_right  Output: reconstructed primitive states at the face
/// @param dx  Cell size (for heaviside evaluation on reconstructed phi)
inline void muscl_reconstruct_primitive_face(
    __local float* l_rho, __local float* l_p, __local float* l_u, __local float* l_v,
    __local float* l_phi,
    int idx_m2, int idx_m1, int idx_C, int idx_p1,
    float dx,
    Conserved* U_left, Conserved* U_right,
    Primitive* P_left, Primitive* P_right
) {
    // --- Left state: extrapolated from center cell (idx_C) toward face ---
    float s_rho = minmod(l_rho[idx_C] - l_rho[idx_m1], l_rho[idx_p1] - l_rho[idx_C]);
    float s_p   = minmod(l_p[idx_C]   - l_p[idx_m1],   l_p[idx_p1]   - l_p[idx_C]);
    float s_u   = minmod(l_u[idx_C]   - l_u[idx_m1],   l_u[idx_p1]   - l_u[idx_C]);
    float s_v   = minmod(l_v[idx_C]   - l_v[idx_m1],   l_v[idx_p1]   - l_v[idx_C]);
    float s_phi = minmod(l_phi[idx_C] - l_phi[idx_m1], l_phi[idx_p1] - l_phi[idx_C]);

    // INTERFACE FLATTENING: Zero slopes when any stencil cell crosses the material
    // interface (phi changes sign). This prevents MUSCL from extrapolating across
    // the 270:1 density / 6:1 sound-speed mismatch, which causes pressure blowup.
    // Detection by sign change is robust even before redistancing makes phi a
    // signed-distance function (initial phi may be ±1 everywhere).
    // Ref [4] Saurel & Abgrall (1999) — flatten reconstruction near interfaces.
    bool xface_L = (l_phi[idx_m2] * l_phi[idx_m1] <= 0.0f) ||
                   (l_phi[idx_m1] * l_phi[idx_C]  <= 0.0f) ||
                   (l_phi[idx_C]  * l_phi[idx_p1] <= 0.0f);
    if (xface_L) { s_rho = 0; s_p = 0; s_u = 0; s_v = 0; s_phi = 0; }

    float rho_L = l_rho[idx_C] + 0.5f * s_rho;
    float p_L   = l_p[idx_C]   + 0.5f * s_p;
    float u_L   = l_u[idx_C]   + 0.5f * s_u;
    float v_L   = l_v[idx_C]   + 0.5f * s_v;
    float phi_L = l_phi[idx_C] + 0.5f * s_phi;

    // --- Right state: extrapolated from neighbor cell (idx_m1) toward face ---
    float t_rho = minmod(l_rho[idx_m1] - l_rho[idx_m2], l_rho[idx_C] - l_rho[idx_m1]);
    float t_p   = minmod(l_p[idx_m1]   - l_p[idx_m2],   l_p[idx_C]   - l_p[idx_m1]);
    float t_u   = minmod(l_u[idx_m1]   - l_u[idx_m2],   l_u[idx_C]   - l_u[idx_m1]);
    float t_v   = minmod(l_v[idx_m1]   - l_v[idx_m2],   l_v[idx_C]   - l_v[idx_m1]);
    float t_phi = minmod(l_phi[idx_m1] - l_phi[idx_m2], l_phi[idx_C] - l_phi[idx_m1]);

    // INTERFACE FLATTENING for right state (same sign-change test)
    bool xface_R = (l_phi[idx_m2] * l_phi[idx_m1] <= 0.0f) ||
                   (l_phi[idx_m1] * l_phi[idx_C]  <= 0.0f) ||
                   (l_phi[idx_C]  * l_phi[idx_p1] <= 0.0f);
    if (xface_R) { t_rho = 0; t_p = 0; t_u = 0; t_v = 0; t_phi = 0; }

    float rho_R = l_rho[idx_m1] - 0.5f * t_rho;
    float p_R   = l_p[idx_m1]   - 0.5f * t_p;
    float u_R   = l_u[idx_m1]   - 0.5f * t_u;
    float v_R   = l_v[idx_m1]   - 0.5f * t_v;
    float phi_R = l_phi[idx_m1] - 0.5f * t_phi;

    // --- Positivity floors on reconstructed primitives ---
    if (rho_L < MIN_RHO) rho_L = MIN_RHO;
    if (rho_R < MIN_RHO) rho_R = MIN_RHO;
    if (p_L < 1.0f) p_L = 1.0f;
    if (p_R < 1.0f) p_R = 1.0f;

    // --- Convert primitives to conserved via inverse EOS ---
    // E = internal_e + KE = (p * denom + pinf_term) + 0.5*rho*(u²+v²)
    // CRITICAL: Use CELL-CENTERED phi (not extrapolated phi_L/phi_R) for alpha.
    // Extrapolating phi across the interface creates inconsistent (rho, alpha)
    // combos — e.g. uranium density with hydrogen alpha → insane pressure.
    // Ref [4] Saurel & Abgrall (1999) — use cell-centered volume fractions.
    float alpha2_L = heaviside(l_phi[idx_C], dx);
    float alpha1_L = 1.0f - alpha2_L;
    float g1_m1 = GAMMA_1 - 1.0f, g2_m1 = GAMMA_2 - 1.0f;
    float denom_L = alpha1_L / g1_m1 + alpha2_L / g2_m1 + EPS;
    float pinf_L  = alpha1_L * GAMMA_1 * PI_1 / g1_m1 + alpha2_L * GAMMA_2 * PI_2 / g2_m1;
    float ie_L    = p_L * denom_L + pinf_L;
    float ke_L    = 0.5f * rho_L * (u_L*u_L + v_L*v_L);

    float alpha2_R = heaviside(l_phi[idx_m1], dx);
    float alpha1_R = 1.0f - alpha2_R;
    float denom_R = alpha1_R / g1_m1 + alpha2_R / g2_m1 + EPS;
    float pinf_R  = alpha1_R * GAMMA_1 * PI_1 / g1_m1 + alpha2_R * GAMMA_2 * PI_2 / g2_m1;
    float ie_R    = p_R * denom_R + pinf_R;
    float ke_R    = 0.5f * rho_R * (u_R*u_R + v_R*v_R);

    // Sound speed from EOS (frozen)
    float gamma_avg_num_L = alpha1_L * GAMMA_1 / g1_m1 + alpha2_L * GAMMA_2 / g2_m1;
    float c_sq_L = (gamma_avg_num_L / denom_L) * (p_L + pinf_L / gamma_avg_num_L) / (rho_L + EPS);
    P_left->c  = sqrt(max(c_sq_L, 1.0f));
    P_left->p  = p_L; P_left->u = u_L; P_left->v = v_L;
    U_left->rho_a1 = rho_L * alpha1_L;
    U_left->rho_a2 = rho_L * alpha2_L;
    U_left->ru     = rho_L * u_L;
    U_left->rv     = rho_L * v_L;
    U_left->E      = ie_L + ke_L;
    U_left->phi    = l_phi[idx_C];  // Cell-centered phi, not extrapolated

    float gamma_avg_num_R = alpha1_R * GAMMA_1 / g1_m1 + alpha2_R * GAMMA_2 / g2_m1;
    float c_sq_R = (gamma_avg_num_R / denom_R) * (p_R + pinf_R / gamma_avg_num_R) / (rho_R + EPS);
    P_right->c = sqrt(max(c_sq_R, 1.0f));
    P_right->p = p_R; P_right->u = u_R; P_right->v = v_R;
    U_right->rho_a1 = rho_R * alpha1_R;
    U_right->rho_a2 = rho_R * alpha2_R;
    U_right->ru     = rho_R * u_R;
    U_right->rv     = rho_R * v_R;
    U_right->E      = ie_R + ke_R;
    U_right->phi    = l_phi[idx_m1];  // Cell-centered phi, not extrapolated
}

/// @brief Main fluid update kernel — 2nd-order MUSCL reconstruction + SSP-RK2 time integration.
/// Each 8x8 tile is processed by 64 threads; 2-cell halo (12x12) loaded cooperatively into local memory.
///
/// SSP-RK2 (Heun's method / Shu-Osher) for dt from t^n to t^{n+1}:
///   Stage 1 (predictor): U* = U^n + dt * L(U^n)           [forward Euler half-step]
///   Stage 2 (corrector): U^{n+1} = 0.5*U^n + 0.5*(U* + dt*L(U*))
/// where L(U) = -(1/dx)*(F_xp - F_xm + F_yp - F_ym) is the spatial operator.
/// Ref: Shu, C.-W. (1988). "Total Variation Diminishing Time Discretizations."
///      SIAM J. Sci. Stat. Comput. 9(6), 1073-1084.
///
/// @param stage 0 = predictor (writes U*), 1 = corrector (reads U^n + U*, writes U^{n+1})
__kernel void update_fluid(
    __global const float4* restrict in_U_vec,
    __global const float2* restrict in_E_vec,
    __global float4* restrict out_U_vec,
    __global float2* restrict out_E_vec,
    __global float*  restrict out_p,
    __global float*  restrict out_c,
    __global const float4* restrict old_U_vec,
    __global const float2* restrict old_E_vec,
    int width,
    int height,
    float dx,
    float dt,
    int b_diagnostics,
    int stage
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
    __local float l_p[LOCAL_SIZE];
    __local float l_u[LOCAL_SIZE];
    __local float l_v[LOCAL_SIZE];
    __local float l_c[LOCAL_SIZE];
    __local float l_rho[LOCAL_SIZE];  // Total density for primitive MUSCL reconstruction

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

    // VACUUM MOMENTUM KILL + COMPUTE PRIMITIVES (cached — avoids 10× redundant EOS per cell)
    for (int t = lx; t < total_items; t += threads) {
        // Kill momentum in vacuum cells BEFORE caching primitives
        if (l_rho_a1[t] + l_rho_a2[t] < MIN_RHO) {
            l_ru[t] = 0.0f;
            l_rv[t] = 0.0f;
        }
        Conserved U;
        U.rho_a1 = l_rho_a1[t];
        U.rho_a2 = l_rho_a2[t];
        U.ru     = l_ru[t];
        U.rv     = l_rv[t];
        U.E      = l_E[t];
        U.phi    = l_phi[t];
        Primitive P = get_primitive(U, dx);
        l_p[t] = P.p;
        l_u[t] = P.u;
        l_v[t] = P.v;
        l_c[t] = P.c;
        l_rho[t] = U.rho_a1 + U.rho_a2;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // COMPUTE — MUSCL 2nd-order reconstruction + Rusanov flux + SSP-RK2
    #define GET_U(idx) (Conserved){l_rho_a1[idx], l_rho_a2[idx], l_ru[idx], l_rv[idx], l_E[idx], l_phi[idx]}
    int pixels_to_process = TILE_W * TILE_H;
    
    for (int k = lx; k < pixels_to_process; k += threads) {
        int ty = k / TILE_W;
        int tx = k % TILE_W;
        int cx = tx + HALO;  // center in local memory (HALO=2 offset)
        int cy = ty + HALO;
        int c_idx = cy * LOCAL_W + cx;

        // --- MUSCL face reconstruction ---
        // For each of 4 faces around cell (cx,cy), reconstruct left/right states
        // using 2-cell stencil: needs indices at cx-2, cx-1, cx, cx+1, cx+2 etc.
        //
        // Face indexing (local memory with HALO=2):
        //   West face  (cx-1/2): left from cell (cx-1), right from cell (cx)
        //   East face  (cx+1/2): left from cell (cx),   right from cell (cx+1)
        //   South face (cy-1/2): left from cell (cy-1), right from cell (cy)
        //   North face (cy+1/2): left from cell (cy),   right from cell (cy+1)

        Conserved U_L_west, U_R_west, U_L_east, U_R_east;
        Conserved U_L_south, U_R_south, U_L_north, U_R_north;
        Primitive P_L_west, P_R_west, P_L_east, P_R_east;
        Primitive P_L_south, P_R_south, P_L_north, P_R_north;

        // West face: between (cx-1) and (cx). Left from (cx-1), right from (cx).
        muscl_reconstruct_primitive_face(l_rho, l_p, l_u, l_v, l_phi,
            cy*LOCAL_W + (cx-2), cy*LOCAL_W + (cx-1), cy*LOCAL_W + cx, cy*LOCAL_W + (cx+1),
            dx, &U_L_west, &U_R_west, &P_L_west, &P_R_west);

        // East face: between (cx) and (cx+1). Left from (cx), right from (cx+1).
        muscl_reconstruct_primitive_face(l_rho, l_p, l_u, l_v, l_phi,
            cy*LOCAL_W + (cx-1), cy*LOCAL_W + cx, cy*LOCAL_W + (cx+1), cy*LOCAL_W + (cx+2),
            dx, &U_L_east, &U_R_east, &P_L_east, &P_R_east);

        // South face: between (cy-1) and (cy). Left from (cy-1), right from (cy).
        muscl_reconstruct_primitive_face(l_rho, l_p, l_u, l_v, l_phi,
            (cy-2)*LOCAL_W + cx, (cy-1)*LOCAL_W + cx, cy*LOCAL_W + cx, (cy+1)*LOCAL_W + cx,
            dx, &U_L_south, &U_R_south, &P_L_south, &P_R_south);

        // North face: between (cy) and (cy+1). Left from (cy), right from (cy+1).
        muscl_reconstruct_primitive_face(l_rho, l_p, l_u, l_v, l_phi,
            (cy-1)*LOCAL_W + cx, cy*LOCAL_W + cx, (cy+1)*LOCAL_W + cx, (cy+2)*LOCAL_W + cx,
            dx, &U_L_north, &U_R_north, &P_L_north, &P_R_north);

        // Rusanov flux at each face using MUSCL-reconstructed states
        Conserved F_xm = compute_flux_from_primitive(U_L_west,  U_R_west,  P_L_west,  P_R_west,  dx, dt, 1.0f, 0.0f);
        Conserved F_xp = compute_flux_from_primitive(U_L_east,  U_R_east,  P_L_east,  P_R_east,  dx, dt, 1.0f, 0.0f);
        Conserved F_ym = compute_flux_from_primitive(U_L_south, U_R_south, P_L_south, P_R_south, dx, dt, 0.0f, 1.0f);
        Conserved F_yp = compute_flux_from_primitive(U_L_north, U_R_north, P_L_north, P_R_north, dx, dt, 0.0f, 1.0f);

        float dtdx = dt / dx;

        // Spatial operator L(U) = -(1/dx)*(dF_x + dF_y)
        Conserved dU;
        dU.rho_a1 = -dtdx * (F_xp.rho_a1 - F_xm.rho_a1 + F_yp.rho_a1 - F_ym.rho_a1);
        dU.rho_a2 = -dtdx * (F_xp.rho_a2 - F_xm.rho_a2 + F_yp.rho_a2 - F_ym.rho_a2);
        dU.ru     = -dtdx * (F_xp.ru     - F_xm.ru     + F_yp.ru     - F_ym.ru);
        dU.rv     = -dtdx * (F_xp.rv     - F_xm.rv     + F_yp.rv     - F_ym.rv);
        dU.E      = -dtdx * (F_xp.E      - F_xm.E      + F_yp.E      - F_ym.E);

        // Non-conservative phi advection (uses center-cell velocity from cache)
        Primitive Pc;
        Pc.p = l_p[c_idx]; Pc.u = l_u[c_idx]; Pc.v = l_v[c_idx]; Pc.c = l_c[c_idx];
        Conserved Uc = GET_U(c_idx);
        float grad_phi_x = (l_phi[cy*LOCAL_W + (cx+1)] - l_phi[cy*LOCAL_W + (cx-1)]) / (2.0f * dx);
        float grad_phi_y = (l_phi[(cy+1)*LOCAL_W + cx] - l_phi[(cy-1)*LOCAL_W + cx]) / (2.0f * dx);
        dU.phi = -dt * (Pc.u * grad_phi_x + Pc.v * grad_phi_y);

        // SSP-RK2 time integration
        // Stage 0 (predictor): U* = U^n + dt * L(U^n) = U^n + dU
        // Stage 1 (corrector): U^{n+1} = 0.5*U^n + 0.5*(U* + dt*L(U*))
        //                       = 0.5*U^n + 0.5*(in_U + dU)  where in_U = U* from stage 0
        Conserved U_new;
        if (stage == 0) {
            // Predictor: U* = U^n + dU
            U_new.rho_a1 = Uc.rho_a1 + dU.rho_a1;
            U_new.rho_a2 = Uc.rho_a2 + dU.rho_a2;
            U_new.ru     = Uc.ru     + dU.ru;
            U_new.rv     = Uc.rv     + dU.rv;
            U_new.E      = Uc.E      + dU.E;
            U_new.phi    = Uc.phi    + dU.phi;
        } else {
            // Corrector: U^{n+1} = 0.5*U^n + 0.5*(U* + dU*)
            // old_U = U^n (saved before predictor), in_U = U* (predictor output)
            int gx_out = tile_origin_x + tx;
            int gy_out = tile_origin_y + ty;
            int g_idx = gy_out * width + gx_out;
            float4 old_U = old_U_vec[g_idx];
            float2 old_E = old_E_vec[g_idx];
            U_new.rho_a1 = 0.5f * old_U.x + 0.5f * (Uc.rho_a1 + dU.rho_a1);
            U_new.rho_a2 = 0.5f * old_U.y + 0.5f * (Uc.rho_a2 + dU.rho_a2);
            U_new.ru     = 0.5f * old_U.z + 0.5f * (Uc.ru     + dU.ru);
            U_new.rv     = 0.5f * old_U.w + 0.5f * (Uc.rv     + dU.rv);
            U_new.E      = 0.5f * old_E.x + 0.5f * (Uc.E      + dU.E);
            U_new.phi    = 0.5f * old_E.y + 0.5f * (Uc.phi    + dU.phi);
        }
        
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
            if (b_diagnostics) {
                Primitive P_new = get_primitive(U_new, dx);
                out_p[out_idx] = P_new.p;
                out_c[out_idx] = P_new.c;
            }
        }
    }
}

/// @brief Redistance phi to signed-distance property via Godunov-type HJ scheme — preserves interface sharpness.
/// Uses upwind gradient selection based on sign of original phi; iterated externally for ~5 steps.
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
