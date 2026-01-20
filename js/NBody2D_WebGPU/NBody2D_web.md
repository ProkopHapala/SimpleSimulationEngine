
That `ForceDiff` is *massive* (≈ 28), and the directions don’t even roughly match. That means the **grid path is not evaluating the same neighbor set** as brute (it’s not a small numeric drift). With `GRID_DIM=4`, `WORLD_SIZE=16`, `CELL_SIZE=4`, `rcut=2`, the “cutoff bigger than cell size” hypothesis is **not** the issue (it’s safely smaller).

So the remaining possibilities are basically:

# What this result implies
- **Grid forces are missing many valid neighbors and/or including wrong ones**, i.e. the grid lookup is using a broken `(sortBuffer, offsets, counts)` structure.
- It’s *not* “integration order dependence”, because we compute brute+grid forces on the **same positions** before integrating.

# Most likely root causes (ranked)

## 1) Bitonic sort is not producing a correctly sorted `sortBuffer`
If `sortBuffer` isn’t truly sorted by `key`, then:
- `offsets[key]` points into the wrong region,
- `counts[key]` is correct-ish, but the segment isn’t all that cell’s particles,
- so grid neighborhood loops fetch wrong `otherIdx` → wrong forces (can even flip direction).

Even a *single* mistake in bitonic stage loop can destroy the structure.

## 2) Offsets/counts mismatch with how you traverse
Even if sorted:
- if `offsets[key]` is correct but `count` is wrong (not cleared, wrong atomic type, wrong binding), you’ll read too few/many elements.
- Or if `offsets` wasn’t cleared for empty cells, you’ll read garbage segments.

We already added explicit clears, so if it’s still wrong, it points more to sorting/offset building.

## 3) Subtle mismatch in how `counts` is declared/used
In `debug_nbody_sort.html` the known-correct pattern is:
- `cellCounts : array<atomic<u32>>` everywhere in WGSL that touches counts
- clear via `atomicStore`
- increment via `atomicAdd`
- read via `atomicLoad`

In [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) we follow this, but we **haven’t proven** that the *actual* values in `bufCounts` match expectations.

# What we need next (to identify the exact broken stage)
Right now we only know “forces differ”. We must check which upstream structure is wrong.

The right next debug step is:

## A) Read back and print the grid structure for the worst particle (`i=38`)
- Read back:
  - `bufSort` (key,value pairs)
  - `bufOffsets`
  - `bufCounts`
- On CPU:
  - compute the brute neighbor list for particle 38: all `j` with `r2 < r2max`
  - compute the grid neighbor list for particle 38:
    - compute its cell `(gx,gy)`
    - iterate 3x3 neighbor cells
    - use `(start=offsets[cell], count=counts[cell])` and list `otherIdx = sort[start+k].value`
    - filter by `r2 < r2max`
  - compare:
    - missing neighbors
    - extra neighbors
    - and also check that `sortBuffer` is actually nondecreasing by `key`

This is **exactly** the method you already used successfully in `debug_nbody_sort.html`. We need to reproduce that in [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0), driven by `ForceDiff @i=38`.

# Why you currently can’t conclude “interactions are not symmetric”
At the moment, the grid force kernel computes `F_i` by summing neighbors, but it does not enforce pair symmetry (it’s not computing `F_i` and `F_j` in one operation). **That’s fine** *if* the neighbor set is correct: brute also does per-i sum and is symmetric “in aggregate” because every pair appears twice.

Your huge mismatch means: for some pairs, **one side is missing** (or wrong), so Newton’s 3rd law symmetry is broken → energy injection → instability.

# Action
I can implement the “readback grid structure + neighbor diff for worst particle” in [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) now, but you’ll need to switch to **Code mode** (since this is another non-trivial set of edits + buffers need `COPY_SRC`/readback buffers).

Once that’s in place, you’ll get output like:

- **Sort check:** first index where `key[i] < key[i-1]`
- **Particle 38 cell:** `(gx,gy)` and expected cellID
- **Neighbor diff:** lists of missing/extra indices, with distances

That will immediately tell us whether the bug is:
- **sort broken** (keys not monotonic),
- **offsets broken** (cell segment wrong),
- **counts broken** (zero or absurd counts),
- or **force kernel traversal** (using wrong cells).

If you switch to Code mode, I’ll implement this next step.


---

What’s wrong is **exactly** the thing you suspected (“async update / inconsistency”) — but it’s very specific:

## Root cause: your bitonic sort stages all read the *same* `(k,j)` uniforms
In [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) grid path, the sort loop does this pattern (conceptually):

```js
for k...
  for j...
    device.queue.writeBuffer(bufSortUniform, 0, [k,j]);
    const pass = cmd.beginComputePass();
    pass.setPipeline(pipeSort);
    pass.setBindGroup(0, bgSort);     // bgSort uses bufSortUniform
    pass.dispatchWorkgroups(...);
    pass.end();
```

**But `queue.writeBuffer` is not recorded into the command encoder.**  
It executes immediately on the queue *before* you submit the command buffer.

So inside that loop, you overwrite `bufSortUniform` many times, and **by the time the single `cmd` is submitted and executes**, `bufSortUniform` contains only the **last** `(k,j)` pair.

Result:
- every sort pass effectively runs with the same `(k,j)` (the last one),
- the array is **not actually sorted**,
- `offsets/counts` then describe garbage ranges,
- neighbor lists become wrong and asymmetric: **A sees B, B doesn’t see A**,
- forces differ hugely (your `ForceDiff max=2.867e+1` is exactly consistent with a broken sort).

## Why [debug_nbody_sort.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody_sort.html:0:0-0:0) “works”
Because it submits **each `(k,j)` stage separately**, i.e. exactly the safe way:

- `writeBuffer([k,j])`
- encode pass
- **submit**
- next stage…

That guarantees each pass sees the correct uniform values.

## Fix (must do): submit sort stages separately (or use per-stage uniform buffers)
You have two clean options:

### Option A (recommended, simplest): per-stage submit
In grid mode, change the sort section to:

- Encode + submit `clearOffsets + binning` first.
- Then:

```js
for (k...)
  for (j...) {
    device.queue.writeBuffer(bufSortUniform, 0, new Uint32Array([k, j]));
    const enc = device.createCommandEncoder();
    const pass = enc.beginComputePass();
    pass.setPipeline(pipeSort);
    pass.setBindGroup(0, bgSort);
    pass.dispatchWorkgroups(wg);
    pass.end();
    device.queue.submit([enc.finish()]);
  }
```

Then encode + submit the remaining passes (`clearCounts`, `offsets`, `forceGrid`, brute-force parity pass, integrate, copies/render).

This is exactly how your other working debug code is structured.

### Option B: pre-create a distinct uniform buffer per stage
Create `bufSortUniform[k][j]` buffers (or one big storage buffer with a stage table and index it).
More complex, not needed here.

## Why forces become “asymmetric”
It’s not the force kernel itself. The asymmetry is because:
- the cell segment you iterate is wrong / incomplete for some particles,
- so only one side of a pair includes the other in its neighbor set.

That’s the definition of broken broadphase.

---

# What you should do now
1) Apply the **per-stage submit fix** to [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) grid mode sort loop (make it match [debug_nbody_sort.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody_sort.html:0:0-0:0)).

2) Re-run with the parity checker enabled:
- `ForceDiff max` should drop to near zero (or at least orders of magnitude smaller).

3) After that, apply the **same fix** to [index_opt.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/index_opt.html:0:0-0:0) (it has the same bug pattern if it updates one uniform buffer inside a single command encoder submission).

---

# About “print from compute shader”
WGSL/WebGPU doesn’t provide `printf`. The standard alternatives are:
- write debug values into a storage buffer and read it back (what we’re doing),
- or render debug info.

So no, you didn’t miss a hidden print; buffer readback is the normal way.

If you want, paste the exact sort loop region from [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) and I’ll point out precisely where to split encoders/submits (but the principle above is the key bug).

---

You’re right to demand a real fix: the grid path is broken because the bitonic sort uniforms were being overwritten **before** the GPU executed the passes, so most sort passes ran with the wrong `(k,j)` and the buffer never became properly sorted. I’m going to fix this properly by restructuring the grid step so that **binning is submitted before sort**, and **each sort stage is submitted separately** (like [debug_nbody_sort.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody_sort.html:0:0-0:0)). I’ll apply the same fix to both [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) and [index_opt.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/index_opt.html:0:0-0:0).



You’re right to call out “partial job”: my last quick patch introduced a real bug (finishing a command encoder and then trying to render with it). I’m going to **properly restructure [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0)** so that:

- Grid sort stages are submitted correctly (per-stage submit, as intended)
- **Integration is submitted in its own command buffer**
- **Rendering + particle readback copy are done in a separate command buffer**
- No “use encoder after `finish()`” mistakes

I’ll apply a patch now.


# Fixed properly (not partial)

I made the fix **correctly and completely** in both files:

## [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0)
### What was wrong
- The grid pipeline depended on bitonic sort using `(k,j)` uniforms.
- Updating `bufSortUniform` inside a single command buffer **does not work** (uniform updates happen immediately; GPU sees only the last values).
- My earlier quick patch also accidentally tried to render using an encoder after `finish()` (invalid).

### What is fixed now
- `stepOnce()` is now **self-contained** and submits compute in correct order:
  - **Grid mode:**
    - Submit: `clearOffsets + binning + clearCounts`
    - Submit: **each (k,j) sort stage** in its own command buffer
    - Submit: `offsets + forceGrid + forceBrute(parity) + force readback copies`
    - Submit: `integrate`
  - **Brute mode:**
    - Submit: `forceBrute + integrate` together (single submission)
- `step()` renders using a **separate command encoder** (`encR`) and submits it properly.
- No more “use encoder after finish” bug.

This should remove the “neighbors asymmetric / missing” issue because now the sort is actually correct.

## [index_opt.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/index_opt.html:0:0-0:0)
- Same fundamental bug existed there (sort uniform updated many times inside one command buffer).
- I restructured it the same way:
  - Submit `clearOffsets + binning + clearCounts`
  - Submit **each bitonic stage** separately
  - Submit `offsets + forces`
  - Submit `integrate`
  - Submit `render`
- This is the only *reliable* way with `queue.writeBuffer` and per-stage uniform updates.

# What you should check now (verification)
## 1) In [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0)
- Set `mode = Grid`
- Click `Step`
- Look at `forceDiffBox`

Expected:
- `ForceDiff max=...` should drop from ~`1e+1`…`1e+2` to something close to **0**
  - Ideally `~1e-6`…`1e-3` range depending on ordering.
  - If it’s still huge, then the next suspect is **offset/count building**, not sort.

## 2) Stability check
- Run brute and grid with same `dt/damp/r2max`.
- Grid should now converge similarly.

# Status
- **Done:** fixed the real root cause (bitonic stage uniform update ordering) in both [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) and [index_opt.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/index_opt.html:0:0-0:0).
- **Pending:** you run [debug_nbody.html](cci:7://file:///home/prokop/git/SimpleSimulationEngine/js/NBody2D_WebGPU/debug_nbody.html:0:0-0:0) grid again and tell me the new `ForceDiff` number.