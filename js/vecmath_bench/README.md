# Vector Math Benchmark

This directory contains benchmarks for the high-performance `Vec3` library.

## Files

- `benchmark_vec3.js`: Node.js benchmark script comparing Object-based vs Array-based Vec3 implementations.
- `index.html`: WebGL visualization of the N-body simulation using `Vec3.js`.

## Running the Benchmark (Node.js)

To run the performance benchmark in your terminal:

```bash
node js/vecmath_bench/benchmark_vec3.js
```

## Running the Visualization (Browser)

1. Ensure you are running a local web server (e.g., `python3 -m http.server` or the project's `run_web_server.sh`).
2. Open `http://localhost:8000/js/common_js/vecmath_bench/index.html` (adjust port if necessary).
3. You should see a particle simulation.
4. Open the browser console (F12) to see performance logs.

## Benchmark Results (Node.js)

**Configuration:** N=1000 bodies, 100 steps, O(N^2) interactions.

| Implementation | Mode | Time (ms) | Steps/sec | Interactions/sec |
| :--- | :--- | :--- | :--- | :--- |
| **Vec3Obj** | **Manual** | **~270 ms** | **~371** | **3.71e+8** |
| Vec3Obj | Methods | ~410 ms  | ~244 | 2.44e+8 |
| Vec3Lst | Manual  | ~431 ms  | ~232 | 2.32e+8 |
| Vec3Lst | Methods | ~522 ms  | ~192 | 1.92e+8 |
| Vec3Arr | Manual  | ~302 ms  | ~331 | 3.31e+8 |
| Vec3Arr | Methods | ~1306 ms | ~77  | 7.65e+7 |
| THREE   | Manual  | ~271 ms  | ~369 | 3.69e+8 |
| THREE   | Methods | ~403 ms  | ~248 | 2.48e+8 |
| GLM     | Manual  | ~451 ms  | ~222 | 2.22e+8 |
| GLM     | Methods | ~746 ms  | ~134 | 1.34e+8 |


**Conclusions:**
1. **Objects (`Vec3.js` and `THREE.Vector3`) are the fastest**. V8 optimizes objects with stable shapes (`x, y, z`) extremely well.
2. **Manual optimization pays off**: Unrolling vector operations manually yields ~30% performance gain over method calls.
3. **TypedArrays have overhead**: `gl-matrix` (which uses TypedArrays) and our `Vec3_arr` are slower than plain objects for individual vector math.
4. **THREE.js is excellent**: Its `Vector3` implementation is almost identical in performance to our custom `Vec3Obj`. Our library has the advantage of being lightweight and having specialized fused operations (`addMul`, `setLincomb`) which THREE.js might lack or require multiple calls to achieve.
