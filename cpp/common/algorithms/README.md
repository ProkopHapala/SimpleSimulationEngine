# algorithms

Header-only algorithm library for simulation and mesh processing. No external dependencies beyond `macroUtils.h` and standard C++.

- **arrayAlgs.h** — Hand-written array algorithms for small-to-medium arrays (n=10..10000): prune (two-pointer compaction by predicate), insertSort (in-place and permutation variants, preferred for n<32 or nearly-sorted), quickSort (Lomuto partition with permutation), objects2cells (CSR counting-sort for Buckets), binSearchBetween (interval lookup for interpolation), exportNonZero (dense-to-sparse extraction). No heap allocation in hot paths.
- **binarySwitch.h** — Bitwise-dispatched level lookup in O(log n). Builds result index one bit at a time from MSB to LSB by testing `levels[i | bit]`. Fixed iteration count makes it branch-predictable and GPU-suitable. Template version (binarySwitchT) unrolls at compile time. bakeSwitchTable pre-computes dense lookup from sparse mapping.
- **sweep.h** — Sweep-and-prune broadphase collision detection. Sorts objects by interval xmin, walks forward with early break when xmax exceeded — O(N·k) for sparse scenes. collideSelf for single set, collideCross for two-set (two-pointer walk, both directions). Span is 8 bytes (float xmin/xmax) for cache efficiency.
- **dataprocess1D.h** — 1D signal processing: bisecNoise1D (fractal midpoint displacement noise), runningMax/runningMin (windowed extrema), findMax. Used for terrain generation and data analysis.
- **synchronize.h** — Double-buffer synchronizer template. Ptr2<T> holds paired pointers (read/write), Synchronizer copies y→x on read and x→y on write. For ping-pong buffer management in time-stepped simulations. (Note: currently incomplete/prototype.)
- **main.cpp** — Test driver for arrayAlgs. Sorts a sample array using quickSort with permutation, prints before/after. Run with: `g++ main.cpp -o test && ./test`
