# Solve Particle Collisions and Spring constraints in 2D

The goal is to make simple demonstration in python how to simulate collisions between molecules composed of spherical atoms.
For simplicity we do it in 2D as it is easier to draw and debug.

## Equality and unequality constraints

- The bonds between atoms are equality constraints, characterized by equlibrium bond length $r_0$ and spring constant $k$. 
   - force $f = k(r-r_0)$ 
- The collisions between non-bonded atoms and the wall are unequality constraints, which operate only when $r>r_0$.
   - force is also  $f = -k(r-r_0)$ but only when $r>r_0$.
   - when $r<r_0$ the force is zero ( which can be formally expressed as k=0).

## Position based approach

- Naieve simulation of dynamics can be done e.g. by Leapfrog integrator by following updates:
   ```
   f = get_force(p);
   v = v + f * dt;
   p = p + v * dt;  
   ```
- However, such approach is not numerically stable especially when stiffness of the spring $k$ is high. This require small time step $dt$, and hence slow (computationally expensive) simulation.
- Therefore it is better idea to use Position-Based dynamiscs (PBD) or its newer variants such as projected dynamics (PD).
- In this approach instead of evaluating the force we evaluate the predicted position of the particle, resp. prediced displacement of particles to minimize the violation of the constraints. 
    - The violation of the constraints is actually the force, i.e. displacement scaled by the stiffness. However, thanks to linearity we know how far to displace the particles.
- Hoewver to ensure physically realistic dynamics we need to upadate velocity after positions were update to satisfy the constraints.

## Position based dynamics

1. Predicted position of particle $i$ is 
   $v' = v + f dt$
   $p' = p + v' dt$
2. solve for new positions $p''$ which satisfy constrains ( i.e. minimize strain) 
    $K p'' = f(p')$
3. update velocity $v$ using the new corrected positions (satisfying the constraints) 
   $v = (p''-p)/dt$
   $p=p''$

## Averaging the displacements

- if particle is affected by single constraint, then the displacement is simply taken to satisfy the constraint.
- if particle is affected by multiple constraints, then we need some weighted linear combination of these displacements $d_i = \sum_j c_{ii} d_{ij}$.
- How to calculate the weights $c_{ij}$?
  - Naturally the more stiff constrains should have more weight.
  - Intuitivelly we can write 
  $$
  d_i = \frac{ \sum_j k_{ij} d_{ij} }{ \sum_j k_{ik} } 
  $$
  therefore
  $$
  c_{ij} = \frac{k_{ij}}{ \sum_j k_{ik} } 
  $$

### From Projective-Dynamics 

Projective dynamics puts this stiffnes weighted displacement on regorous grounds and further introduce inertia into the equation by adding position $p'$ predicted by dynamics update and the scalled by $M/dt^2$. 

Matrix equation of projective dynamics can be written as system of N linear equationd. Each equation is centered aroung point $p_i$. 
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

This can be acctually written as 
$$
p_i = \frac{ \sum_j ( K_{ij} (  d_{ij} - p_j )) }{ \sum_j K_{ij} }
$$

resp. if we introduce back the intertia

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

## Collisions offset

We indroduce the offset $D$ to the collision detection, which includes the collision into the constrain solver only if $r<r_0-D$. This has two reasons:
- Collisions can be solved as linear springs well if we are deep in repulsion (penetration) stage.
- Also we do not need to solve collisions by constrain solver as long as we are in mild-repulsion stage. In this regime we can safely use non-linar dynamical relaxation based of forces and velocities with reasonably large time step $dt$.

## Update of collision neighbors

- We solve both bond-lenght and non-bonded collisions byt sparse approach using neighbor-list for each atom, without explicitly forming any matrix. We solve the costraints iteratively e.g. by Jacobi or Gauss-Seidel method.
- While bonding topology is fixed, non-bonded collision neighbor list is updated on the fly.
- If we consider constrain solver as a linear matrix equation, we should keep the neighbor list constant for all iteration of the solver.
- However it makes sense to experiment with updating the neighbor in between iterations of the solver, in the hope that it will introduce non-linear effects inherent to collision (inequality constraints).

## Neighborhood search

- To do this collision neighbor list is updated efficiently we should pre-screen and short lit possible neighbors of given particle $i$ before we start the costrain solver loop. 
- In fact, assuming the time step and velocity of particles is small, we can do this costly $O(n^2)$ operation only once per several dynamical steps.  
- Ideally we can share thist short-list of all possible neighbors between multiple nearby particles, because thanks to trinagle inequality, we can be sure that particles $i,j$ which are close enough ($r_{ij}<R_g$) to each other are also close enough $r_{ik}<R_c+R_g$ to particle $k$ if $r_{kj}<R$.
- The best way how to exploit this idea is group the particles into chunks of rize $R_g$ based on their proximity, and list all particles with distance less than $R_c$.
- There are many spatial data structures how to do this.
- For our demonstration in 2D with densly packed molecules two approches seems most suitable
  1. grid-based approach where particles are first assigned by grid-cells by rounding their position, and then searching neihborst in the adjacent cells.
  2. using each molecule as the group, and searhing neighbors in the molecules which centers of mass are closer then $R_g+R_c$.
- Initially we will implement the molecule-based  grouping approach as it is the most simple one, and it has little overhead for large sparse world.