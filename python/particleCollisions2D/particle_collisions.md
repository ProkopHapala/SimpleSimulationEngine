# Solve Particle Collisions and Spring constraints in 2D

The goal is to make a simple demonstration in Python of how to simulate collisions between molecules composed of spherical atoms.
For simplicity, we do it in 2D as it is easier to draw and debug.

## Equality and Inequality Constraints

- The bonds between atoms are equality constraints, characterized by equilibrium bond length $r_0$ and spring constant $k$. 
   - force $f = k(r-r_0)$ 
- The collisions between non-bonded atoms and the wall are inequality constraints, which operate only when $r>r_0$.
   - force is also  $f = -k(r-r_0)$ but only when $r>r_0$.
   - when $r<r_0$ the force is zero (which can be formally expressed as k=0).

## Position-Based Approach

- Naive simulation of dynamics can be done e.g. by Leapfrog integrator with following updates:
   ```
   f = get_force(p);
   v = v + f * dt;
   p = p + v * dt;  
   ```
- However, such approach is not numerically stable especially when stiffness of the spring $k$ is high. This requires small time step $dt$, and hence slow (computationally expensive) simulation.
- Therefore it is better to use Position-Based dynamics (PBD) or its newer variants such as projected dynamics (PD).
- In this approach, instead of evaluating the force we evaluate the predicted position of the particle, resp. predicted displacement of particles to minimize the violation of the constraints. 
    - The violation of the constraints is actually the force, i.e. displacement scaled by the stiffness. However, thanks to linearity we know how far to displace the particles.
- However, to ensure physically realistic dynamics we need to update velocity after positions were updated to satisfy the constraints.

## Position-Based Dynamics

1. Predicted position of particle $i$ is 
   $v' = v + f dt$
   $p' = p + v' dt$
2. solve for new positions $p''$ which satisfy constraints ( i.e. minimize strain) 
    $K p'' = f(p')$
3. update velocity $v$ using the new corrected positions (satisfying the constraints) 
   $v = (p''-p)/dt$
   $p=p''$

## Averaging the Displacements

- If a particle is affected by a single constraint, then the displacement is simply taken to satisfy the constraint.
- If a particle is affected by multiple constraints, then we need some weighted linear combination of these displacements $d_i = \sum_j c_{ij} d_{ij}$.
- How to calculate the weights $c_{ij}$?
  - Naturally the more stiff constraints should have more weight.
  - Intuitively we can write 
  $$
  d_i = \frac{ \sum_j k_{ij} d_{ij} }{ \sum_j k_{ij} } 
  $$
  therefore
  $$
  c_{ij} = \frac{k_{ij}}{ \sum_j k_{ij} } 
  $$

### From Projective-Dynamics 

Projective dynamics puts this stiffness weighted displacement on rigorous grounds and further introduces inertia into the equation by adding position $p'$ predicted by dynamics update and scaled by $M/dt^2$. 

Matrix equation of projective dynamics can be written as system of N linear equations. Each equation is centered around point $p_i$. 
$$
( m_i/dt^2 + \sum_j K_{ij} )  p_i   +    \sum_j (K_{ij} p_j)  = p'_i m_i/dt^2 + \sum_j (K_{ij} d_{ij})  
$$

where:
 - $p_i$ is the position of particle $i$ which we search for
 - $p'_i$ is the predicted position of particle $i$ by dynamics ignoring internal constraint forces
 - $m_i$ is the mass of particle $i$
 - $K_{ij}$ is the spring constant of bond between particles $i$ and $j$ (or the collision stiffness)
 - $d_{ij}$ is the displacement of particle $i$ which would satisfy the constraint between particles $i$ and $j$ ignoring all other constraints.

The above equation can be simplified by expressing the diagonal of the matrix $K_ii$ and right hand side $b_i$ as
 - $ K_{ii} = m_i/dt^2 + \sum_j K_{ij}$
 - $ b_i    = ( m_i/dt^2)p'_i + \sum_j (K_{ij} d_{ij})$

One can notice that if we ignore the dynamical term (considering inertial term $M/dt^2 = 0$ ), then the system of equation can be written as
$$
( \sum_j K_{ij} ) p_i - \sum_j (K_{ij} p_j)  = \sum_j (K_{ij} d_{ij})  
$$

This can be actually written as 
$$
p_i = \frac{ \sum_j ( K_{ij} (  d_{ij} - p_j )) }{ \sum_j K_{ij} }
$$

resp. if we introduce back the inertia

$$
p_i = \frac{ \sum_j ( K_{ij} (  d_{ij} - p_j )) + m_i/dt^2 p'_i }{ \sum_j K_{ij} + m_i/dt^2 }
$$

This absolutely makes sense, even more if we express $p'_{ij} = d_{ij} + p_j$ as the predicted position of particle $i$ which satisfy the constraint with particle $j$. Then we have
$$
p_i = \frac{ \sum_j ( K_{ij} p'_{ij}) + m_i/dt^2 p'_i }{ \sum_j K_{ij} + m_i/dt^2 }
$$ 

Which is basically just weighted average of predicted positions of particles due to dynamics and due to all the constraints.

NOTE: notice that last equation is exactly in form of Jacobi update 

$$
p_i = \frac{ b_i -\sum_j ( K_{ij} p_j ) }{ K_{ii} }
$$ 

## Collisions Offset

We introduce the offset $D$ to the collision detection, which includes the collision into the constraint solver only if $r<r_0-D$. This has two reasons:
- Collisions can be solved as linear springs well if we are deep in repulsion (penetration) stage.
- For mild-repulsion stage, we can safely use non-linear dynamical relaxation based on forces and velocities with reasonably large time step $dt$.

Note: The stiffness $K$ of the constraint solver at the transition point $r = R_c - D$ should match the stiffness of the non-linear inter-atomic non-covalent potential (e.g., Lennard-Jones) used in the dynamics.

## Update of Collision Neighbors

- We solve both bond-length and non-bonded collisions by sparse approach using neighbor-list for each atom, without explicitly forming any matrix. We solve the constraints iteratively e.g. by Jacobi or Gauss-Seidel method.
- While bonding topology is fixed, non-bonded collision neighbor list is updated on the fly.
- If we consider the constraint solver as a linear matrix equation, we should keep the neighbor list constant for all iterations of the solver.
- However, it makes sense to experiment with updating the neighbor in between iterations of the solver, in the hope that it will introduce non-linear effects inherent to collision (inequality constraints).

## Neighborhood Search

- To efficiently update the collision neighbor list, we should pre-screen and short-list possible neighbors of given particle $i$ before starting the constraint solver loop. 
- In fact, assuming the time step and velocity of particles is small, we can do this costly $O(n^2)$ operation only once per several dynamical steps.  
- Ideally, we can share this short-list of all possible neighbors between multiple nearby particles, because thanks to triangle inequality, we can be sure that particles $i,j$ which are close enough ($r_{ij}<R_g$) to each other are also close enough $r_{ik}<R_c+R_g$ to particle $k$ if $r_{kj}<R$.
- The best way how to exploit this idea is to group the particles into chunks of size $R_g$ based on their proximity, and list all particles with distance less than $R_c$.
- There are many spatial data structures how to do this.
- For our demonstration in 2D with densely packed molecules, two approaches seem most suitable:
  1. grid-based approach where particles are first assigned to grid-cells by rounding their position, and then searching neighbors in the adjacent cells.
  2. using each molecule as the group, and searching neighbors in the molecules which centers of mass are closer than $R_g+R_c$.
- Initially, we will implement the molecule-based grouping approach as it is the simplest one, and it has little overhead for large sparse worlds.







# Scratchpad


now the error we should plot is the lengh of residual from jacobi_iteration_sparse()

r[i]     = b[i] +sum_j - Aii[i]*x[i]

this is basically 
r = b - Ax

if we use formulas 

A_{ii} = \sum_j K_{ij}
b_i    = \sum_j (K_{ij} d_{ij})

we can say

y=Ax is
y_i = A_{ii} p_i + \sum_j K_{ij} p_j  

r_i = b_i - y_i =   \sum_j (K_{ij} d_{ij}) - A_{ii} p_i - \sum_j K_{ij} p_j 
r_i =    \sum_j ( K_{ij} d_{ij} )  -  \sum_j K_{ij} p_i - \sum_j K_{ij} p_j   
r_i =    \sum_j K_{ij} ( d_{ij} - p_i + p_j )
r_i =    \sum_j K_{ij} ( p'_ij - p_i   )

where   
p'_ij = d_{ij} + p_j is the predicted position of point p_i due to contrain with point p_j
and
d_{ij} = dir_{ij} * r0 (resp. r_c) is the optimal vectro between the two points

Therefore the residual r_i should reflect how far is each point p_i from its optimal contrained position p'_ij therefore it should reflect the fact that we do not satisfy the constrian betwen the two atoms  