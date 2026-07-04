# Ray Tracing Implementation Comparison

## Basic Idea

Both implementations calculate radiosity coupling between surface elements:
- CPU version uses object-oriented C++ with explicit loops
- OpenCL version parallelizes calculations using GPU kernels

## Pseudocode

### CPU Version
```pseudo
for each surface element i:
    for each surface element j:
        calculate geometric coupling
        if coupling > threshold:
            calculate occlusion
            update coupling matrix
```

### OpenCL Version
```pseudo
kernel makeRadiosityCouplings:
    for each work item (parallel surface element i):
        for each surface element j:
            calculate geometric coupling
            if coupling > threshold:
                calculate occlusion
                update coupling matrix
```

## Key Differences

- OpenCL uses float4 vectors for better GPU utilization
- CPU version has more detailed geometric calculations
- OpenCL version handles memory differently (global buffers)
- CPU version includes more debug visualization options

## Detailed Implementation Comparison

### Triangle Intersection

#### CPU Implementation (TriangleRayTracer.h)
- Exact ray-triangle intersection in [getOcclusion()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/cpp/common/dynamics/TriangleRayTracer.h:171:4-187:5) method (lines 142-160)
- Uses `rayInTriangle()` method to check intersection
- Processes obstacles sequentially in a for loop
- Binary result (1 for intersection, 0 for no intersection)
- No distance calculation or fuzzyness
- Example:
```cpp
for(int i=0; i<triangleObstacles.size(); i++) {
    if(triangleObstacles[i].rayIn(ray0,hX,hY)) {
        return 1.0;
    }
}
```

### Acceleration Structures

#### CPU Implementation

- No explicit spatial acceleration structures
- Processes all obstacles sequentially
- Relies on CPU cache for performance
- Simple for loop through all obstacles
```cpp
for(int i=0; i<triangleObstacles.size(); i++)
```

#### OpenCL Implementation
- Work group processing (lines 898-961)
- Chunked obstacle processing

```cl
for(int i0=0; i0<ns.y; i0+= nL)
```

Local memory caching for better memory access patterns
```cl
LOC[iL*3] = obstacles[i3];
```

- Two-stage processing:
    - Groups of points bounded within spheres
    - Groups of triangles bounded within larger triangles
- Uses barrier synchronization for parallel processing


## Triangle Intersection Comparison

### CPU Implementation (TriangleRayTracer.h)
- Uses exact ray-triangle intersection
- Checks if ray intersects triangle using rayInTriangle method
- Returns binary result (1 for intersection, 0 for no intersection)
- Processes obstacles sequentially
- No distance calculation or fuzzyness

### OpenCL Implementation (radiosity.cl)
- Uses local memory caching for obstacle triangles
- Implements parallel processing using work groups
- Includes surface index checking to skip self-intersections
- Processes obstacles in chunks for better memory access patterns
- Uses barrier synchronization between processing stages

## Acceleration Structures

### CPU Implementation
- No explicit spatial acceleration structures
- Processes all obstacles sequentially
- Relies on CPU cache for performance
- Doesn't group obstacles into larger blocks

### OpenCL Implementation
- Uses work groups to process chunks of obstacles
- Caches obstacle triangles in local memory
- Implements two-stage processing:
  1. Groups of points bounded within spheres
  2. Groups of triangles bounded within larger triangles
- Uses barrier synchronization for parallel processing
- Skips self-intersections using surface indices