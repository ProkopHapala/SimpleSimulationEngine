I’m going to make [RARFF_2.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/RARFF_2.html:0:0-0:0) behave like the **relaxation** mode in the C++ reference (no inertial “jerking”), and also ensure the **grid path can’t silently break like in NBody2D.md**. Concretely: default to **brute**, switch integration to an **overdamped/gradient-descent** step (like `RARFF_SR::move()` / `moveMDdamp` but without accumulating velocities), and generate **non-overlapping initial positions** on a lattice instead of random (random starts cause huge forces and collapse). Then we can re-enable grid once the relaxed dynamics look right.

I’ll patch [RARFF_2.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/RARFF_2.html:0:0-0:0) to make the dynamics *relax* like the C++ “move/gradient descent” mode (no inertial velocities), reduce initial overlap by placing atoms on a lattice, and start in **Brute** mode by default. This directly targets the “jerking/collapsing” you’re seeing.
## Fixed to match your goal: **stable relaxation / bond formation** (no jerking)

I applied the key stabilization changes directly in [js/NBody2D_WebGPU/RARFF_2.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/RARFF_2.html:0:0-0:0):

### 1) Default is now **Brute force**
- The UI `<select>` defaults to **Brute**
- `useGrid` now starts as `false`
This avoids the “NBody2D.md” class of issues (broken neighbor set) unless you explicitly switch to Grid.

### 2) Integrator changed to **overdamped relaxation** (no inertial dynamics)
In `integrateWGSL` I removed the “velocity + damping + symplectic-ish” update and replaced it with a **gradient-descent style** step:

- `pos += force * (dt*damping)`
- quaternion updated directly from **torque** (as a relaxation step)
- `vel` and `angVel` are forced to **zero each step**

This matches the “relax” behavior you want (like C++ `move()`-style), and avoids oscillations/jerks that come from inertia + strong short-range forces.

### 3) Initial particle placement fixed (avoids catastrophic overlaps)
Random initial positions cause huge forces → collapse/jerking.  
Now particles start on a **3D lattice** with spacing `~1.8*rbond0`, then wrapped into the box. Velocities/angVel are still zero.

## About the NBody2D.md problem (are we affected?)
- That document describes a **grid-sort uniform update hazard** (all sort stages seeing the last `(k,j)`), which breaks the neighbor set and injects energy.
- In [RARFF_2.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/RARFF_2.html:0:0-0:0) you are already using the **safe approach**: each bitonic stage is submitted separately (so it does *not* silently reuse the same uniform).
- Also: since we now default to **Brute**, you’re safe by default even if grid regresses.

## What to do now
1) Hard reload [RARFF_2.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/RARFF_2.html:0:0-0:0) (cache-bust `?v=9`)
2) Leave it in **Brute** initially.
3) If you still see “fishy” behavior, tell me which failure mode:
- **[collapse]** everything clumps into one blob
- **[explosion]** atoms fly apart
- **[frozen]** almost no motion / no bonding
- **[still oscillating]** you still see jerky motion despite overdamped mode

Then we tune just 2 things:
- **dt** (smaller for stability)
- **aMorse / bMorse / rcut** (bond strength + range)

### Status
- **Done:** default brute, overdamped relax integrator, lattice initialization (major stability fixes).
- **Next (optional):** make Grid safe+fast again (and/or add 2-atom [testEF](cci:1://file:///home/prokop/git/SimpleSimulationEngine/cpp/sketches_SDL/Molecular/test_RARFF_SR.cpp:74:0-100:1) validation mode + picking).