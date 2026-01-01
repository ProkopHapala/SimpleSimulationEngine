# Pressure Debugging Retrospective (EulerianImpacFluid)

This document summarizes the issues encountered with low/negative pressure in the Uranium projectile and the fixes applied. It serves as a checklist for future work to avoid similar pitfalls.

## Root Causes Identified
1. **Volume Fraction vs. Mass Fraction Mix-up**
   - The 5-equation EOS requires volume fractions (alpha), but the kernel initially derived alpha from partial densities (mass fractions), making Uranium appear as a tiny volume fraction in interface cells. This drove pressure to the floor.
   - Fix: Use `alpha2 = heaviside(phi, dx); alpha1 = 1 - alpha2` in the kernel EOS (OpenCL) and in host diagnostics.

2. **Phi Advection Was Conservative**
   - `phi` was updated with a conservative flux, smearing the signed distance field and de-correlating volume fractions from partial densities.
   - Fix: Switch `phi` update to non-conservative advection: `phi_t + u*phi_x + v*phi_y = 0` with central gradients.

3. **Inconsistent Initialization Between Phi and Partial Densities**
   - Partial densities were set via hard masks (sharp), while energy used a smooth heaviside, causing EOS to see mismatched `rho_mix` vs. `alpha` near the interface.
   - Fix: Initialize `rho_a1 = RHO_H2 * (1 - H(phi))`, `rho_a2 = RHO_U * H(phi)`, and energy from the same `H(phi)` so interface cells are consistent and isobaric at p0.

4. **Missing Momentum/Energy in Projectile Init**
   - Projectile kinetic energy was initially omitted; total energy was too low, collapsing pressure.
   - Fix: Add KE = 0.5 * rho_mix * v^2 into E during init; set momentum ru = rho_tot * v_smooth.

5. **GPU Diagnostics Not Populated at t=0**
   - Pressure/sound buffers read back as zeros before first compute.
   - Fix: Run `solver.step(0.0)` then `solver.queue.finish()` before first frame readback.

6. **Precision Sensitivity with Large PI_2 (40 GPa)**
   - Uranium’s stiffened gas PI makes pressure extremely sensitive to small internal-energy deficits. Float32 rounding and interface inconsistencies led to p floors.
   - Fixes that mitigated this:
     - Consistent alpha from phi (both kernel and host diagnostics).
     - Non-conservative phi advection.
     - Isobaric smooth init using the same H(phi) for rho_a*, E, and velocity.
     - Host-side double precision diagnostics to verify E − KE − pinf is positive.

## Key Code Changes (Highlights)
- **Kernel (OpenCL)**
  - EOS: `alpha2 = heaviside(phi, dx); alpha1 = 1 - alpha2`.
  - Phi update: non-conservative advection (central gradients), not fluxed.
  - Diagnostic outputs: pressure `p` and sound speed `c` stored to GPU buffers.

- **Driver (Python)**
  - Added buffers for `p` and `c`; returning them in `get_data()`.
  - After init: `solver.step(0.0)` + `queue.finish()` to populate diagnostics.

- **Test Script Initialization (Python)**
  - Isobaric, smooth init: `rho_a1 = RHO_H2 * (1 - H(phi))`, `rho_a2 = RHO_U * H(phi)`, `ru = rho_tot * v_smooth`, `E = p0*denom + pinf + KE`, all using the same heaviside.
  - Host diagnostics use the same phi-derived alphas; KE in float64 for checks.

## Symptoms Observed & Resolution
- **Symptom:** Pressure in Uranium near zero/1 Pa while Hydrogen near p0.
- **Resolution:** Align alpha with phi in EOS, smooth/init-consistent rho_a*, add KE, non-conservative phi advection, and sync GPU diagnostics.

- **Symptom:** GPU pressure read as zeros at t=0.
- **Resolution:** Run `step(0.0)` and `queue.finish()` before first readback.

- **Symptom:** NaNs/negative internal energy.
- **Resolution:** Clamping KE vs. total E; internal energy floor; consistent init; diagnostics in float64 to verify positivity.

## Other Problems Encountered (NaNs / Units / Floors)
- **Unit mismatches and CFL:** Early defaults used dt too large (e.g., 4e-7) and mixed unit assumptions; switching to SI-consistent params and tighter dt (1e-8) stabilized CFL-driven blowups.
- **Pressure floor masking root cause:** Raising MIN_P to huge values hid the real EOS inconsistency; proper fix is consistent alpha/phi/rho and KE inclusion, not a large floor.
- **Sound speed / PI_2 sensitivity:** Incorrect or inconsistent sound speed calculations and very large PI_2 (40 GPa) amplified small internal-energy deficits into NaNs/negatives; fixed by consistent EOS and precision checks.
- **Missing KE in projectile init:** Total energy too low (internal only) caused immediate pressure collapse and NaNs; adding KE to E fixed this.
- **Phi/rho mismatch at interface:** Sharp masks for rho with smooth phi caused mixed cells to violate EOS assumptions, leading to NaNs/low pressure; fixed by scaling rho_a* with the same H(phi) used for energy and EOS.
- **Diagnostics lag:** Reading GPU p/c before kernel execution returned zeros, confusing debugging; fixed with a zero-step + queue finish.

## Lessons / Takeaways for Developers
1. **Always derive volume fractions from phi for the EOS.** Mass fractions at interfaces will break mixed-material pressure.
2. **Phi is not a conserved density.** Update phi with advection, not conservative fluxes.
3. **Initialization must be self-consistent.** Use the same heaviside for rho_a*, E, and velocity; include kinetic energy for moving projectiles.
4. **Synchronize GPU diagnostics.** After init, run a zero-step and finish the queue before plotting GPU-side p/c.
5. **Check with higher precision.** When PI is huge (stiffened gas), small float32 errors can collapse pressure; validate with float64 diagnostics.
6. **Interface sensitivity.** Smooth transitions (heaviside epsilon) should match between phi, rho_a*, and E to maintain isobaric starts.

## Current Status
- Uranium pressure now initializes near p0 and remains positive; diagnostics show GPU p in expected ranges (no immediate floor to 1 Pa).
- Remaining work: further reduce long-term pressure drift, and consider CFL/PI tuning for stability at very high Mach.
