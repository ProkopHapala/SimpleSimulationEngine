Here’s the math framing of the polar used in [AeroSurf.h](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/dynamics/AeroSurf.h:0:0-0:0):

**Inputs**
- `û`: unit airflow vector in the surface frame.
- `ca = cos(α) = û · n̂` where `n̂` is local spanwise “up” (grot.c is negative forward, grot.b is up).
- `sa = sin(α) = û · t̂` where `t̂` is local chord/upwind tangent (grot.b).
- Stall shaping uses `|sa|`.

**Stall blending**
- `wS = trashold_cub(|sa|; sStall, sStall + wStall)` (smoothstep from 0 to 1 across the stall band).
- `mS = 1 − wS`.

**Drag coefficient**
\[
CD = CD_0 + \big( m_S \, dCD \, |sa| + w_S \, dCDS \big)\, |sa|
\]

**Lift coefficient (with sign handling)**
- If `ca < 0`, flip: `ca = -ca`, `sa = -sa` (keeps lift pointing consistently with effective AoA).
\[
CL = \big( m_S \, dCL + w_S \, dCLS \, ca \big)\, sa
\]

**Force construction**
- Dynamic pressure factor (without ½ρ): `prefactor = |v_air|^2 * area`.
- Lift direction is airflow-normal component:
\[
\hat{\ell} = \frac{t̂ - sa\,û}{\|t̂ - sa\,û\|}
\]
(implemented as `airUp = grot.b – sa*û`, normalized).
- Total force:
\[
\vec{F} = CL \cdot \hat{\ell} \cdot prefactor \;+\; CD \cdot û \cdot prefactor
\]
Applied at the panel’s world-space position; torque comes from `craft->apply_force`.

**Parameters**
- `CD0` base drag; `dCD` linear drag slope pre-stall; `dCDS` post-stall drag slope.
- `dCL` lift slope pre-stall; `dCLS` post-stall lift slope (reduced).
- Stall onset `sStall`, width `wStall`.
- Cutoff: forces skipped if `|v_air|^2 < lowSpeedCuoff`.

This captures how lift/drag are blended from attached flow to stalled, based on `sin`/`cos` of incidence between airflow and the surface normal/tangent.