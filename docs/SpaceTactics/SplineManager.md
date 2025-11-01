# `SplineManager.h`

This header file provides the `SplineManager` class, a tool for creating, managing, and evaluating multi-dimensional Hermite splines. These splines are used extensively in `SpaceTactics` to define the smooth trajectories and thrust profiles of ships.

## `SplineManager` Class

The `SplineManager` class handles a set of control points and their associated time values to define one or more smooth curves.

### Public Members

*   `int n`: The number of control points.
*   `int m`: The number of independent dimensions (or curves) managed by the spline.
*   `double* ts`: An array of time values corresponding to each control point.
*   `double** CPs`: A 2D array storing the values of the control points for each dimension.
*   `double** dCPs`: An optional 2D array storing explicit derivative values at each control point. If not provided, derivatives are inferred from the control points.

### Public Methods

*   `void allocate(int m_, int n_, bool derivs)`: Allocates memory for the spline data, including `n_` control points and `m_` dimensions. If `derivs` is true, memory is also allocated for explicit derivatives.
*   `void deallocate(int n_, int m_)`: Frees all memory allocated for the spline data.
*   `double evalIt(int ip, double t, ..., double* val, double* dval, double* ddval)`: Evaluates the spline at a given time `t`, which is known to be within the time interval `[ts[ip], ts[ip+1]]`. This is the core evaluation function.
    *   It can compute the value (`val`), the first derivative (`dval`), and the second derivative (`ddval`) of the spline.
    *   It uses the Hermite basis functions to ensure a smooth curve that passes through the control points.
*   `int eval(double t, int i0, ..., double* val, double* dval, double* ddval)`: A convenience wrapper that first performs a binary search to find the correct time interval for `t` (starting the search from index `i0`) and then calls `evalIt` to perform the evaluation.
*   `int evalUniform(double t0, double t1, int n, double* val, double* dval, double* ddval)`: Evaluates the spline at a series of uniformly spaced time steps between a start and end time.
*   `int removePoint(int i)`**: Removes the control point at index `i` from the spline.
*   `int insertPoint(double t, int i)`**: (Placeholder) Intended to insert a new control point into the spline at a specific time `t`.
