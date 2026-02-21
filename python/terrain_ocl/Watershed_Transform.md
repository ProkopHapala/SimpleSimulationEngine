

## USER

OK more or less works now.  But perhaps we need to clearly define what we actally want. The choice of gates and sinks is a bit arbitrary. In fact what we want to do is something else.

We should go downhill, every pixel say which pixel is lower (outflow), the sink is any pixel which does not have any lower nighbor (outflow),
the id of outflow pixel propagate back up the hill, until it reach some ridge (saddle) where boundaries of two sink ids meet, that is basin boundary. This define boundary between basins coresponding to each sink. 

then another kernel goes over all pairs of neighbor pixels and find the minimal height at which we have boundary of each sing has minimum (not sure what is best way, maybe propagate again down the hill takeing minimum from neighboring pixels?)

this will generate perfect segmentation of the map into basins and pathways between basins. We can implement it at fist in global memory ignoring the tile boundaries and local memory, since that complicates algorithm and hinder the physics/math. Only when we debug it in the simple global version we will try to reintroduce tiles and local memory. 

think it through before we start to write the kernels.

---

## Gemini-3-Pro

This is the correct pivot. You are describing the **Watershed Transform** (specifically the "Drop Analysis" or "Steepest Descent" method).

This approach is topologically superior because:
1.  **Physics-based:** We don't force "one sink per 16x16 tile". Sinks appear where the terrain actually dictates.
2.  **Natural Boundaries:** Basins are defined by ridges, not grid lines.
3.  **Path Definition:** A path between two basins is simply: *Start at the Saddle Point (Gate) -> Flow down to Sink A, then Flow down to Sink B.*

Here is the design for the **Global Memory / Global Logic** implementation. We will ignore tiles for a moment to ensure the math produces perfect rivers and ridges.

### The Algorithm

1.  **Kernel 1: Flow Calculation (Steepest Descent)**
    For every pixel, look at the 8 neighbors. Point to the one that produces the steepest drop (taking distance $\sqrt{2}$ into account).
    *   If no neighbor is lower, the pixel is a **Sink**.
    *   *Result:* A `flow_map` where every pixel points to its downstream parent.

2.  **Kernel 2: Basin Labeling (Pointer Jumping)**
    We need to know which Sink every pixel belongs to.
    *   Init: `basin_id = (is_sink ? my_index : -1)`
    *   Iterate: `basin_id[i] = basin_id[ flow_map[i] ]`
    *   Repeat until convergence. This effectively collapses the entire tree so every pixel knows its "Ultimate Drain".
    *   *Result:* `basin_id_map`.

3.  **Kernel 3: Gate Discovery (Ridge Detection)**
    We run a kernel over edges (horizontal and vertical pairs).
    *   If `basin_id[x,y] != basin_id[x+1,y]`: We found a **Ridge**.
    *   The "Pass Height" is `max(height[x,y], height[x+1,y])`.
    *   We store this candidate gate.
    *   *Reduction:* On the CPU (or a reduction kernel), we group these candidates by the pair `(basin_A, basin_B)` and keep only the one with the **Minimum Pass Height**.

---

#### 1. `terrain_global.cl`
#### 2. `verify_global.py`


### Why this is better:
1.  **No more zig-zags:** The flow kernel calculates the vector of steepest descent. If the terrain is smooth, the path is smooth.
2.  **No more arbitrary tiles:** Sinks are where they physically exist.
3.  **Perfect Connectivity:** By definition, every basin is connected to its neighbors via the lowest point on their shared ridge.

Shall we proceed with implementing this **Global Watershed** approach?