# TrussDynamics_d.h

This file defines the `TrussDynamics_d` class, which implements a dynamic simulation of a truss structure. It includes functionalities for setting up the truss, applying forces, solving the equations of motion using various numerical methods (including Projective Dynamics and Conjugate Gradient), handling collisions, and providing picking/selection capabilities. It also incorporates features for linearized truss simulation and FIRE optimization.

---

## Includes

- `datatypes.h`: Defines fundamental data types used throughout the project.
- `Vec2.h`: Defines the `Vec2` class for 2D vector operations.
- `Vec3.h`: Defines the `Vec3` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3x3 matrix operations.
- `quaternion.h`: Defines the `Quat4d` and `Quat4f` classes for quaternion operations, used for representing rotations.
- `VecN.h`: Defines the `VecN` class for N-dimensional vector operations.
- `Buckets.h`: Defines the `Buckets` class for spatial partitioning, used for collision detection.
- `raytrace.h`: Defines functions for ray tracing, used for picking and collision detection.
- `geom3D.h`: Defines geometric primitives and functions for 3D geometry.
- `Interfaces.h`: Defines interfaces for interaction with external systems or solvers.
- `CG.h`: Defines the `CGsolver` class for solving linear systems using the Conjugate Gradient method.
- `arrayAlgs.h`: Defines various array-based algorithms.
- `SparseMatrix.h`: Defines the `SparseMatrix` class for storing and manipulating sparse matrices.
- `SparseMatrix2.h`: Defines the `SparseMatrix2` class, an alternative implementation for sparse matrices tailored for specific needs within the truss dynamics simulation, such as Cholesky decomposition.
- <stdio.h>: Standard input/output library.
- <string.h>: Standard string manipulation library.
- `testUtils.h`: Includes utilities for testing and debugging the code.

---

## Free functions

### Force Calculation
- `springForce( double l, double& f, Quat4d par )`: Calculates the spring force magnitude and potential energy based on length and stiffness.
- `springForce( Vec3d d, Quat4d par )`: Calculates the spring force vector and potential energy based on displacement and stiffness.

### Debugging and Validation
- `checkDist(int n, const Vec3d* vec, const Vec3d* ref, int verb=1, double tol=1e-12 )`: Checks the maximum distance between two sets of 3D vectors.
- `print_vector( int n, double * a, int pitch, int j0, int j1 )`: Prints a section of a vector with a specified pitch.

### Bounding Box Fitting
- `fitAABB( Quat8d& bb, int n, int* c2o, Quat4d * ps )`: Fits an Axis-Aligned Bounding Box to a set of points.
- `fitAABB_edge( Quat8d& bb, int n, int* c2o, int2* edges, Quat4d * ps )`: Fits an Axis-Aligned Bounding Box to a set of edges.

### Bounding Box Updating
- `updatePointBBs(const Buckets& buckets, Quat8d* BBs, Quat4d* points, bool bInit=true)`: Updates bounding boxes based on point positions.
- `updateEdgeBBs(const Buckets& buckets, Quat8d* BBs, int2* edges, Quat4d* points, bool bInit=true)`: Updates bounding boxes based on edge positions.

---
## Types (classes and structs)
---

### class `SmartMixer`

The `SmartMixer` class is a helper class designed to control the blending of iterative solutions using a momentum-like mixing strategy. It allows for a smooth transition between different solution states, which can improve convergence and stability in iterative solvers.

#### properties

##### Blending Parameters
- `b_end`: `float` - The blending factor at the end of the mixing process.
- `b_start`: `float` - Unused.
- `istart`: `int` - Iteration number at which blending starts.
- `iend`: `int` - Iteration number at which blending ends.  Currently Unused.

#### methods

##### Blending Factor
- `get_bmix(int itr)`: Returns the blending factor based on the current iteration number.

---

### class `EdgeVertBond`

The `EdgeVertBond` structure defines the relationship between a vertex and an edge in the truss structure, it can be used to constrain dynamics, or eval forces (e.g. collision).

#### properties

##### Geometric and Constraint Parameters
- `verts`: `Vec3i` - Indices of the vertices forming the edge and the vertex being constrained (edge.a,edge.b,ivert).
- `c`: `double` - Interpolation parameter (0 to 1) that determines the position of the constraint point along the edge.
- `K`: `double` - Stiffness constant used to enforce the constraint between the vertex and the edge.
- `f`: `Vec3d` - The force applied by the constraint.

---

### class `TrussDynamics_d`

The `TrussDynamics_d` class is the core class for simulating the dynamics of a truss structure. It manages the points, forces, velocities, and connections within the truss, and it provides methods for applying forces, solving the equations of motion, handling collisions, and enabling user interaction.

**Inheritance**

- Picker: Inherits picking functionality from the `Picker` class, allowing users to select points or bonds within the simulation.

#### properties

##### Simulation Parameters
- `time`: `double` - The current simulation time.
- `Cdrag`: `double` - Drag coefficient applied to the points.
- `dt`: `double` - The simulation time step.
- `kGlobal`: `double` - Global stiffness constant for the truss structure.
- `damping`: `double` - Damping coefficient applied to the simulation.
- `nSolverIters`: `int` - The number of iterations used in the linear solver.
- `kLinRegularize`: `double` - Force constant for linear regularization.
- `cg_tol`: `double` - The tolerance for the Conjugate Gradient solver.
- `time_LinSolver`: `double` - Time taken by the linear solver.
- `time_cg_dot`: `double` - Time taken by the dot product operation in the Conjugate Gradient solver.
- `linSolveMethod`: `int` - An integer selecting the linear solver method to be used.
- `bApplyResudualForce`: `bool` - Flag indicating whether to apply a residual force after the linear solve.
- `residualForceFactor`: `double` - Factor for scaling the residual force.
- `mass`: `double` - Total mass of the truss structure.
- `F_residual`: `double` - Residual force in the system.
- `maxAcc`: `double` - Maximum acceleration allowed for a point.
- `collision_damping`: `double` - Damping factor applied during collisions.
- `lastNeg`: `int` - Counter for consecutive negative values (used in FIRE algorithm).
- `minLastNeg`: `int` - Minimum number of consecutive negative values before increasing the time step (used in FIRE algorithm).
- `finc`: `double` - Factor for increasing the time step (used in FIRE algorithm).
- `fdec`: `double` - Factor for decreasing the time step (used in FIRE algorithm).
- `falpha`: `double` - Factor for damping reduction (used in FIRE algorithm).
- `dt_max`: `double` - Maximum allowed time step.
- `dt_min`: `double` - Minimum allowed time step.
- `damp_max`: `double` - Maximum allowed damping.
- `ff_safety`: `double` - Safety factor to avoid division by zero.
- `cv`: `double` - A variable related to Conjugate Gradient solver
- `cf`: `double` - A variable related to Conjugate Gradient solver

##### Geometric Properties
- `pos0`: `Vec3d` - Position of the rotating frame's origin.
- `ax`: `Vec3d` - Rotation axis of the rotating frame.
- `cog`: `Vec3d` - Center of gravity of the truss structure.
- `vcog`: `Vec3d` - Velocity of the center of gravity of the truss structure.
- `I`: `Mat3d` - Moment of inertia tensor of the truss structure.
- `L`: `Vec3d` - Angular momentum of the truss structure.
- `torq`: `Vec3d` - Torque acting on the truss structure.
- `hit_pos`: `Vec3d` - Position of a ray hit during picking or collision detection.
- `hit_normal`: `Vec3d` - Normal vector at the point of a ray hit.

##### Point and Force Data
- `nPoint`: `int` - The number of points in the truss structure.
- `points`: `Quat4d*` - An array of `Quat4d` objects storing the position (xyz) and mass (w) of each point.
- `forces`: `Quat4d*` - An array of `Quat4d` objects storing the force (xyz) and energy (w) acting on each point.
- `vel`: `Quat4d*` - An array of `Quat4d` objects storing the velocity of each point.
- `vel0`: `Quat4d*` - An array of `Quat4d` objects storing the old velocity of each point.
- `kFix`: `double*` - An array of force constants for fixed points.
- `ps_cor`: `Vec3d*` - An array of `Vec3d` objects storing the corrected positions of the points after the solver step.
- `ps_pred`: `Vec3d*` - An array of `Vec3d` objects storing the predicted positions of the points before the solver step.
- `ps_0`: `Vec3d*` - An array of `Vec3d` objects storing the initial positions of the points.

##### Bond Data
- `nBonds`: `int` - The number of bonds (edges) in the truss structure.
- `bonds`: `int2*` - An array storing the indices of the two points connected by each bond (edge).
- `strain`: `double*` - An array storing the strain of each bond (edge).
- `maxStrain`: `Vec2d*` - An array storing maximum strain (tensile and compressive) for each bond.
- `kDirs`: `Vec3d*` - Normalized direction of the stick (d/|d|) where d = p1-p0.
- `bparams`: `Quat4d*` - An array of `Quat4d` objects storing bond parameters (l0, kPress, kTens, damp).

##### Neighbor Data
- `nNeighMax`: `int` - The maximum number of neighbors for each point.
- `nNeighTot`: `int` - Total number of neighbor entries (nPoint * nNeighMax).
- `neighs`: `int*` - An array storing the indices of neighboring points for each point.
- `neighBs`: `int2*` - An array storing the indices of neighbor bonds for each point (bond index and neighbor index).
- `neighB2s`: `int*` - An array storing the indices of neighboring bonds for each point.
- `params`: `Quat4d*` - An array of `Quat4d` objects storing neighbor parameters (l0, kP, kT, damp).

##### Linear Solver Data
- `bvec`: `Quat4d*` - Right-hand side vector (b) for the linear system Ap=b in Projective Dynamics, storing internal forces.
- `bvec0`: `Vec3d *` - Backup right-hand side vector for Projective Dynamics.
- `extern_b`: `Quat4f*` - External right-hand side vector for external linear solver.
- `extern_x`: `Quat4f*` - External solution vector for external linear solver.
- `PDmat`: `double*` - The matrix for Projective Dynamics, typically a dense matrix.
- `LDLT_L`: `double*` - The L (lower triangular) matrix from the LDLT decomposition of the PDmat.
- `LDLT_D`: `double*` - The diagonal matrix from the LDLT decomposition of the PDmat.
- `neighsLDLT`: `int*` - Neighbor indices for the LDLT decomposition.
- `nNeighMaxLDLT`: `int` - The maximum number of neighbors considered in the LDLT decomposition.
- `Lsparse`: `SparseMatrix2<double>` - Sparse matrix representation of the L matrix.
- `LsparseT`: `SparseMatrix2<double>` - Sparse matrix representation of the transpose of the L matrix.
- `PDsparse`: `SparseMatrix<double>` - Sparse matrix representation of the PDmat.
- `linsolve_b`: `Vec3d*` - An array of `Vec3d` objects used as the right-hand side vector in the linear solver.
- `linsolve_yy`: `Vec3d*` - An array of `Vec3d` objects used as an intermediate vector in the linear solver.

##### Damping
- `damped_bonds`: `std::vector<int>` - A vector of indices of bonds that are subject to additional damping.
- `damped_points`: `std::unordered_set<int>` - A set of indices of points that are subject to additional damping.

##### Collision Detection
- `nBBs`: `int` - Number of bounding boxes used for collision detection.
- `BBs`: `Quat8d*` - An array of bounding boxes (can be AABB, cylinder, or capsule).
- `pointBBs`: `Buckets` - Buckets for collision detection between points.
- `edgeBBs`: `Buckets` - Buckets for collision detection between edges.
- `faceBBs`: `Buckets` - Buckets for collision detection between faces.
- `pointChunks`: `Buckets` - Chunks for parallelization, these points are copied to local memory when solving one edgeBBs chunk of bonds
- `nFaces`: `int` - Number of faces in the truss structure.
- `faces`: `int4*` - An array storing the indices of points forming triangles or quads, used for ray-tracing and collision detection.
- `nEdgeVertBonds`: `int` - Number of edge-vertex bonds.
- `edgeVertBonds`: `EdgeVertBond*` - An array of `EdgeVertBond` objects representing the connections between vertices and edges.

##### External Solver
- `extern_solve`: `void (*extern_solve)()` - Function pointer to an external linear solver that can be used to solve the system Ap=b.

##### Distortion
- `hbs`: `Quat4f*` - Array storing normalized direction of the stick, and initial distortion of the length from neutral length
- `dpos`: `Quat4f*` - Array storing the distortion of point from neutral postion
- `fdpos`: `Quat4f*` - Array storing the force on distortion of point from neutral postion
- `vdpos`: `Quat4f*` - Array storing the velocity of distortion of point from neutral postion
- `points_bak`: `Quat4d*` - Backup of points for position based dynamics.
- `dls`: `float*` - Distortion of point from neutral postion.

##### Helper Classes
- `mixer`: `SmartMixer` - An instance of the `SmartMixer` class used for blending between different solution strategies in iterative solvers, such as momentum mixing.
- `cgSolver`: `CGsolver` - An instance of the `CGsolver` class used to solve the linear system using the Conjugate Gradient method.

##### User Interaction
- `user_update`: `void (*user_update)(double dt)` - Function pointer to a user-defined update function that is called between iterations.

##### Constant Forces
- `accel`: `Quat4d` - Constant acceleration vector applied to all points.
- `rot0`: `Quat4d` - Center of rotation for applying rotational forces.
- `omega`: `Quat4d` - Angular velocity for the rotating frame.

#### methods

##### Simulation Control
- `set_time_step`: Sets the simulation time step.
- `run`: Runs the simulation for a specified number of iterations.
- `setOpt`: Sets the time step and damping values for the simulation.
- `run_omp`: Runs the simulation with OpenMP for parallelization.

##### Force Calculation
- `getPointForce`: Calculates the force acting on a point due to drag and acceleration.
- `applyForceRotatingFrame_i`: Applies Coriolis and centrifugal forces to a point in a rotating frame.
- `applyForceCentrifug_i`: Applies centrifugal force to a point.
- `applyForceRotatingFrame`: Applies Coriolis and centrifugal forces to the entire truss in a rotating frame.
- `applyForceCentrifug`: Applies Centrifugal forces to the entire truss.
- `evalTrussForce`: Evaluates the truss force.
- `evalTrussForce_neighs`: Evaluates the truss force for a given point using neighbor data.
- `evalTrussForce_neighs2`: Evaluates the truss force for a given point using neighbor data from precomputed neighbors.
- `evalTrussForces_neighs2`: Evaluate truss forces for all points using neighbors2 scheme.
- `evalTrussForces_neighs`: Evaluate truss forces for all points using neighbors scheme.
- `evalTrussForces_bonds`: Evaluates the truss forces based on the bonds (edges).
- `evalBondTension`: Evaluates the tension in each bond.

##### Linear Solver
- `prepare_LinearSystem`: Prepares the linear system for solving, including allocation and matrix assembly.
- `run_LinSolve`: Runs the linear solver for a specified number of iterations.
- `run_Cholesky_omp_simd`: Runs the Cholesky solver with OpenMP and SIMD.
- `run_Cholesky_omp`: Runs the Cholesky solver with OpenMP.
- `solveLinearizedConjugateGradient`: Solve the linearized truss system using conjugate gradient method.
- `dotPD`: Calculates the matrix-vector product Ap.

##### Linearized Truss Simulation
- `prepareLinearizedTruss`: Prepares the linearized truss for simulation.
- `evalTrussForcesLinearized`: Evaluates the forces in the linearized truss.
- `evalTrussForceLinearized_neighs2`: Evaluates the linearized truss force for the point iG.
- `evalTrussForcesLinearized_neighs2`: Evaluates the linearized truss forces for all the points.
- `prepareLinearizedTruss_ling`: Sets up the linearized truss by computing initial forces and stiffness values.
- `dot_Linearized_bonds`: Perform matrix multiplication for linearized bonds.
- `dot_Linearized_neighs2`: Perform matrix multiplication for linearized truss system using neighbor indices.
- `findGuassSeidelPivotingPriority`: Find Guass Seidel pivoting priority.

##### Position Based Dynamics
- `move_dpos`: Move points using forces, and backup the point position.
- `apply_dpos`: Apply displacement on points.
- `update_velocity`: Update the point velocities using the point displacements.
- `constr_jacobi_neighs2_absolute`: Solve the system by constrain distances.
- `constr_jacobi_neighs2_diff`: Constrain dynamics Jacobi iterative difference solver.
- `run_constr_dynamics`: Runs the constrained dynamics simulation.

##### Neighbor Management
- `recalloc`: Reallocates memory for the simulation data structures.
- `reallocFixed`: Reallocates the `kFix` array to match the current number of points.
- `printNeighs`: Prints the neighbor information for a given point.
- `printAllNeighs`: Prints the neighbor information for all points.

##### Collision Detection
- `recallocBBs`: Reallocates memory for bounding box data structures.
- `edgesToBBs`: Assigns edges to bounding boxes based on proximity to points.
- `printBBs`: Prints the bounding box information.
- `evalTrussCollisionImpulses_bonds`: Evaluates the collision impulses between points connected by bonds.

##### Picking
- `getPickedObject`: Returns a pointer to the picked object based on the `mask` value.
- `pick_point_brute`: Performs brute-force picking of a point closest to a ray.
- `pick_bond_brute`: Performs brute-force picking of a bond closest to a ray.
- `pick_BBox`: Picks a bounding box intersected by a ray.
- `pick_nearest`: Picks the nearest object (point or bond) to a ray.
- `pick_all`: Picks all objects (points or bonds) within a certain radius of a ray.

##### Projective Dynamics
- `make_PD_Matrix`: Constructs the Projective Dynamics matrix.
- `make_PDmat_sparse`: Constructs the sparse Projective Dynamics matrix.
- `rhs_ProjectiveDynamics`: Calculates the right-hand side vector for the Projective Dynamics linear system.
- `rhs_ProjectiveDynamics_i`: Calculates the right-hand side vector for Projective Dynamics for a single node `i`.
- `rhs_ProjectiveDynamics_`: Calculates the right-hand side vector for the Projective Dynamics linear system using `rhs_ProjectiveDynamics_i`.
- `updatePD_RHS`: Updates the right-hand side vector for Projective Dynamics.
- `updatePD_dRHS`: Updates right hand side (RHS) of linear system with position differences.

##### Iterative Solvers
- `updateJacobi_lin`: Performs a Jacobi iteration for solving the linear system.
- `updateGaussSeidel_lin`: Performs a Gauss-Seidel iteration for solving the linear system.
- `updateJacobi_fly`: Performs a Jacobi iteration "on-the-fly," calculating the force and diagonal terms within the iteration.
- `updateGaussSeidel_fly`: Performs a Gauss-Seidel iteration "on-the-fly," calculating the force and diagonal terms within the iteration.
- `updateIterativeMomentum`: Update positions using iterative momentum solver.
- `updateIterativeJacobi`: Updates positions using iterative Jacobi method.
- `updateIterativeJacobiDiff`: Solve the displacement from initial position using Jacobi iteration.
- `updateIterativeMomentumDiff`: Solve the displacement from initial position using Momentum iteration.
- `updateIterativeExternDiff`: Solve the displacement from initial position using external linear solver.

##### FIRE Optimization
- `FIRE_update`: Update variable of FIRE algorithm.

##### Helper Functions
- `getLinearBondStiffness`: Returns the linear bond stiffness for a given bond.
- `norm_butFixed`: Calculates the norm of a vector of positions, excluding fixed points.
- `updateInveriants`: Updates the invariants of the truss structure (mass, center of gravity, moment of inertia, angular momentum).
- `cleanForce`: Resets all forces to zero.
- `cleanVel`: Resets all velocities to a specified value.

##### Movement
- `move_GD`: Moves the points using gradient descent.
- `move_i_MD`: Moves the point using Molecular Dynamics.
- `move_MD`: Moves all the points using Molecular Dynamics.

##### Angular Velocity
- `addAngularVelocity`: Adds angular velocity to the points around a specified axis.
- `addAngularVelocity2`: Adds angular velocity to the points around a specified axis, scaled by the distance from axis.