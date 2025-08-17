# Sub‑Pixel Contour Fitting with Compact Bases

This module (`python/subPixelContour/`) fits sharp 2D contours with a small grid of basis functions, enabling efficient, low‑resolution representations of objects with crisp boundaries.

- Use cases: collision masks on grids, binary materials/regions in multi‑component simulations, fonts and glyphs, map features (lakes, walls), implicit shapes for CAD/SDF.
- Goals: accurate zero‑level set (contour) with few grid nodes; robust to irregular sampling; fast evaluation; tunable locality and smoothness.

Main driver: `subPixelContour.py`.


## Quickstart
- Run as a module (requires repository root on PYTHONPATH):
```
export PYTHONPATH=.
python -u -m python.subPixelContour.subPixelContour --opt dyn --niter 400 --plot_opt_log
```
- Or run the script directly:
```
python -u python/subPixelContour/subPixelContour.py --opt dyn --niter 400 --plot_opt_log
```
Dependencies: `numpy`, `matplotlib`.


## Problem Statement
Given sample points `x_p ∈ R^2` and binary reference values `b_p ∈ {−1, +1}` describing the interior/exterior of a shape, we fit an implicit function `F(x)` whose zero‑contour approximates the shape. We construct `F` from compact local bases placed on a coarse `nx×ny` grid of nodes.

Let `B_i(x)` be a basis function centered at grid node `i`. We consider:
- Linear model: `F(x) = Σ_i c_i B_i(x)`.
- Rational (normalized) model: `F(x) = N(x)/D(x)` with `N(x)=Σ_i c_i w_i B_i(x)`, `D(x)=Σ_i w_i B_i(x)`, and non‑negative weights `w_i ≥ 0` (either fixed uniform or learned non‑uniform).

We do not aim to fit the full function away from the contour. Instead, we define a clamped error that focuses learning near `F≈0` (details below).


## Basis Families (kind)
Basis sampling matrix `A ∈ R^{P×N}` has `A[p,i] = B_i(x_p)` at sample `p` and node `i`. Built by `build_A(points, basis_points, kind)` in `subPixelContour.py`. Supported kinds:

- bilinear: `B(dx,dy) = max(0,1−|dx|) · max(0,1−|dy|)` (support in a single grid cell, separable, C0)
- bilinear_power (aka bilinear_p, bilinear_pow): uses `B^p` with exponent `p>0` and target transform (see below)
- bilinear_nu (aka bilinear_nonuniform, bilinear_w): same bilinear `A` but used with the rational model and optional non‑uniform weights
- bspline (cubic B‑spline): separable tensor product of 1D cubic B‑splines with support `|t|<2`
  - 1D cubic B‑spline: for `a=|t|`:
    - if `a<1`: `B(t) = 2/3 − a^2 + 1/2 a^3`
    - if `1≤a<2`: `B(t) = (1/6) (2−a)^3`
- rbf2x2 (radial over 2×2 cells) on normalized radius `r = sqrt(((x−x_i)/h)^2 + ((y−y_i)/h)^2)`. Select kernel via basis kind:
  - `rbf2x2` → r2_linear: `max(0, 1 − r^2)`
  - `rbf2x2_sq` → r2_quadratic: `max(0, 1 − r^2)^2`
  - `wendland` → Wendland C2: `(1 − r)^4_+ (4 r + 1)`

All bases are compact and local; with regular grids, only the 4 corners of the containing cell contribute for 2×2‑support families.


## Target Transform for bilinear_power
When using `kind ∈ {bilinear_power, bilinear_p, bilinear_pow}`, we fit in a transformed space to balance sensitivity:
- Basis is replaced by `B^p` (element‑wise power, `p>0`).
- Targets are transformed before fitting: `b_fit = sign(b) · |b|^p`.
- For visualization only, predictions are inverted: `F_vis = sign(F) · |F|^{1/p}`.

This helps shape the transition near the contour while keeping the local support of bilinear bases.


## Error Model and Clamping (focus on the contour)
We want the zero‑level set `F=0` to match the reference contour, not to fit values far from it. We therefore define a clamped residual per sample:
- `φ(F) = clip(F, −1, 1)`
- `e = φ(F) − b` with `b ∈ {−1, +1}`

The derivative of the clamp is:
- `dφ/dF = 1` if `F ∈ (−1, 1)`, else `0`

Thus only samples with predictions inside `|F|<1` contribute gradients. This focuses learning on a band around the zero‑contour. Implementation details:
- Linear model force/grad: `make_force_linear()` computes `g = A^T (dφ ∘ e)`; the descent force is `f = −g`.
- Clamping is applied consistently in `grad_rational_coeffs()` and `grad_rational_weights()` for the rational model as well.


## Linear Model: Gradient
For `F = A c` with clamped residuals as above:
- `e = φ(F) − b`
- `dφ = 1` on `(−1,1)` and `0` outside
- Gradient: `g = ∂J/∂c = A^T (dφ ∘ e)` for the least‑squares cost `J = 1/2 Σ_p e_p^2`

Code: `make_force_linear()` uses this gradient to provide a model‑agnostic force to relaxators.


## Rational Model: Normalized Combination and Non‑Uniform Weights
To stabilize blending and reduce ghosting, we use a normalized (rational) combination
```
F = N/D,
N = A (w ∘ c),
D = A w,
```
with `w_i ≥ 0` and optional learning of `w` (non‑uniform weights). This ensures local normalization by the sum of contributions.

### Gradient wrt coefficients c
Using clamped error `e` and `dφ` as above, one gets
```
∂F_p/∂c_k = (w_k A[p,k]) / D_p
Let t = (e ∘ dφ) / D   (element‑wise division)
Then g_c = ∂J/∂c = (A^T t) ∘ w
```
Implementation: `grad_rational_coeffs(A,c,b,w, clamp=True)`.

### Gradient wrt weights w (non‑uniform)
The derivative of the normalized model wrt `w_k` is
```
∂F_p/∂w_k = A[p,k] (c_k − F_p) / D_p
Let t = (e ∘ dφ) / D
Then g_w = (A^T t) ∘ c − A^T (t ∘ F)
```
Implementation: `grad_rational_weights(A,c,b,w, clamp=True)`.

We enforce positivity/robustness by clipping learned weights: `w ∈ [wclip_min, wclip_max]` (see `make_force_rational()` and `project()` therein). Learning `w` adapts local normalization, sharpening/softening transitions where needed while preserving locality.


## Model‑Agnostic Relaxators and Logging
Optimization is intentionally generic and separate from model details:
- `gradient_descent_relax(force_fn, ...)` — simple explicit descent in DOF space
- `dynamical_relaxation(force_fn, ...)` — damped dynamics with velocity and opposition checks

Both accept `force_fn(u) -> (f, E)` from model factories and log iterations as a list of tuples `(E, ||f||)` for easy plotting by `plot_opt_log()`.


## List of functions (in `subPixelContour.py`)

### Model setup and basis construction
- `generate_basis_points()` — integer‑grid node positions `(n,2)` used as basis centers.
- `generate_grid_points()` — regular sampling grid and stacked point list for evaluation/plotting.
- `build_A()` — constructs `A[p,i]=B_i(x_p)` for chosen `kind`; handles bilinear power and dispatches to evaluators; non‑uniformity handled later by solvers.
- `evaluate_bilinear_bases()` — separable bilinear `max(0,1−|dx|)·max(0,1−|dy|)` for all points vs nodes.
- `evaluate_bspline_bases()` — 2D cubic B‑spline via tensor product of `cubic_bspline_1d()` in x/y.
- `evaluate_rbf2x2_bases()` — compact radial kernels over a 2×2 neighborhood (`r2_linear`, `r2_quadratic`, `wendland_c2`).

### Reference preparation
- `reference_shape()` — binary circle (+1 inside, −1 outside) used in demos.
- `points_in_any_circle()` — OR of circles; `points_inside_convex_polygon()` — inside test for convex polygon (CCW).
- `reference_from_shapes()` — combine circles OR outside polygon; returns boolean then mapped to ±1.
- `evaluate_reference()` — selects analytic reference or composed masks; returns ±1 values at points.

### Model evaluation and gradients
- `predict_rational()` — normalized prediction `F=N/D` with `N=A@(w∘c)`, `D=A@w`; optional clamp for diagnostics; returns `(F,N,D)`.
- `grad_rational_coeffs()` — gradient of `1/2‖φ(F)−b‖²` wrt `c`: `g=(A^T t)∘w`, with `t=(e·dφ)/D`.
- `grad_rational_weights()` — gradient wrt `w`: `g_w=(A^T t)∘c − A^T(t∘F)` (same clamped band).

### Model‑specific force factories
- `make_force_linear()` — returns `(force, project, u_dim)` for linear clamped model (`u=c`); `force(u)->(−∂J/∂c, ‖e‖)`.
- `make_force_rational()` — returns `(force, project, dim, unpack)` for rational model; packs/unpacks `u=[c,w]` when training `w`, applies weight clipping in `project()`.

### Model‑agnostic optimizers (relaxators)
- `gradient_descent_relax()` — explicit DOF stepper; supports mass scaling/projection; logs `(E,‖f‖)`; early‑stops by `‖f‖`.
- `dynamical_relaxation()` — damped dynamics with velocity and opposition reset; supports mass/projection; logs `(E,‖f‖)`.

### Plotting and diagnostics
- `plot_results()` — reference, model (clamped/full), and error fields; prints RMSE; overlays basis points.
- `plot_opt_log()` — plots optimization traces of `|Error|` and `|Force|` from list‑of‑tuples log (supports legacy dicts).
- `plot_single_basis_2d()` — visualize one basis column of `A` as a 2D image over a grid.
- `plot_single_basis_1d()` — visualize one basis along a 1D line segment through the domain.

### Legacy specialized solvers
- `solve_with_gradient_descent()` / `solve_with_damped_dynamics()` — linear model (pre‑relaxator path).
- `solve_with_gradient_descent_rational()` / `solve_with_damped_dynamics_rational()` — rational model (pre‑relaxator path).


## Minimal API usage (programmatic)
```python
import numpy as np
from python.subPixelContour import subPixelContour as SPC

# grid and samples
nx, ny, res = 10, 10, 120
basis_pts = SPC.generate_basis_points(nx, ny)
xmin, xmax, ymin, ymax = -0.5, nx-0.5, -0.5, ny-0.5
_, _, pts = SPC.generate_grid_points(xmin, xmax, ymin, ymax, res=res)

# reference (circles ∪ outside polygon) mapped to ±1
circles = [(nx*0.35, ny*0.35, 1.8), (nx*0.65, ny*0.65, 1.2)]
poly    = np.array([[0,0],[nx-2,1],[nx-2,ny-2],[1,ny-2]], float)
b = SPC.evaluate_reference(pts, circles=circles, polygon=poly)

# basis and model choice
kind = 'wendland'                 # or 'bilinear', 'bspline', 'rbf2x2', 'rbf2x2_sq'
A = SPC.build_A(pts, basis_pts, kind=kind)

# force factory (rational with non‑uniform weights)
force, project, u_dim, unpack = SPC.make_force_rational(A, b, clamp=True, wclip_min=1e-8, wclip_max=None, train_w=True)
u0 = np.concatenate([np.zeros(A.shape[1]), np.ones(A.shape[1])])

# relaxator
log = []
u = SPC.dynamical_relaxation(force, u0, dt=0.1, damping=0.1, n_iterations=400, ftol=None, log=log, project=project,
                             mass=np.concatenate([np.full(A.shape[1],1.0), np.full(A.shape[1],1.0)]) )
c, w = unpack(u)

# visualization on the same grid
F = SPC.predict_rational(A, c, w=w, clamp=False)[0].reshape(res, res)
ref = SPC.evaluate_reference(pts, circles=circles, polygon=poly).reshape(res, res)
SPC.plot_results(ref, F, [xmin,xmax,ymin,ymax], basis_pts)
SPC.plot_opt_log(log)
```

## CLI Usage (examples)
Run the demo (bilinear, damped dynamics, plot logs):
```
python -u -m python.subPixelContour.subPixelContour --opt dyn --niter 400 --plot_opt_log
```
Try rational with non‑uniform weights and Wendland basis:
```
python -u -m python.subPixelContour.subPixelContour --basis wendland --rational --nonuniform --opt dyn --niter 600 --plot_opt_log
```
Use bilinear power with exponent p (and transformed targets):
```
python -u -m python.subPixelContour.subPixelContour --basis bilinear_power --p 2.0 --opt gd --niter 400 --plot_opt_log
```
Common flags:
- `--size S` or `--nx N --ny M` set grid resolution
- `--res R` sampling resolution per axis
- `--rational`, `--nonuniform` enable normalized model and learning of `w`
- `--dt` (or legacy `--lr`) time step; `--damping` (or legacy `--mu`) for dynamics
- `--ftol`, `--wftol` early stop by force norms
- `--wclip_min`, `--wclip_max` bounds for learned weights
- `--test_const` diagnostic for rational constancy between equal neighbors


## Diagnostics and Advanced Flags
- __`--opt {gd,dyn}`__: choose relaxator (`gd`=gradient descent, `dyn`=damped dynamics).
- __Mass scaling__: `--mass_c`, `--mass_w` set relative update inertia for coefficients vs weights (rational NU). Larger mass = slower updates.
- __Early stop__: `--ftol`, `--wftol` thresholds on force norms for c and w.
- __Constancy diagnostic__: `--test_const` checks F along a line for two equal neighbors (use with `--basis rbf2x2|rbf2x2_sq|wendland`).
- __Legacy params__: `--lr`, `--mu`, `--wlr`, `--wmu` kept for compatibility. Generic relaxators use `--dt` and `--damping`; prefer `--mass_w` to slow weights rather than `--wlr`.


## Practical Notes
- Clamping focuses learning to a band `|F|<1`; outside this band gradients are zero. Choose sampling resolution and band width consistently with your target scale.
- The normalized model mitigates amplitude bias from uneven support coverage; learning `w` (non‑uniform) can adapt local normalization to complex geometries while staying compact.
- Basis choice controls locality/smoothness: bilinear (sharpest), B‑spline (smoother), Wendland (radial, C2).


## What’s New Here
- Compact, local bases with a contour‑centric clamped error that ignores regions far from the zero‑set.
- Rational normalization to stabilize blending, with optional non‑uniform weights learned from data.
- Clean separation of model forces from generic relaxators, enabling easy experimentation.

