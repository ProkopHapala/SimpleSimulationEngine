## Splitting Leap-Frog propagation with constraint solver

Propagation of Newtons equations of motion using Leap-Frog integration scheme is done by following two steps:

1. update velocity at half step: $v_{k+/2} = v_{k-/2} + (f^{tot}_k/m)    dt$
2. update position at full step: $p_{k+1}    = p_k    + v_{k+/2} dt$

However, in position-based projective dynamics we do work with explicit total force $f^{tot} = f^{ext} + f^{con}$, because the internal forces originating from hard constrains $f^{con}$  (e.g. springs, collisions) may be too high which would make the dynamics unstable and necessitate very small time steps.

Instead we propagate dynamical equations with only soft external forces $f^{ext}$ and correct the position by constraint solver. Then we need to update velocity to satisfy the equations of motion with updated position. The whole propagation scheme is done by following two steps:
1. predict velocity at half step: $v'_{k+/2} = v_{k-/2} + (f^{ext}_k/m)   dt$
2. predict position at full step: $p'_{k+1}  = p_k      +  v_{k+/2} dt$
3. correct positions to satisfy constraints: $ p_{k+1} = A^{-1} b( p'_{k+1} )$,
   where the right hand side $b$ is computed from predicted positions $p'_{k+1}$
4. correct velocity at half step: $v_{k+/2} = (p_{k+1} - p_k) / dt$

## Conservation of momentum

Now, the last step $v_{k+/2} = (p_{k+1} - p_k) / dt$ is questionable. We should check if this velocity satisfies the equations of motion and does not violate some laws of physics, in particular the conservation of energy and linear and angular momentum.

We start by discretization of Newton's equations of motion:
$$ 
\frac{d^2}{dt^2} p(t) = f^{tot}(t)/m
$$
By replacing derivatives by finite differences we get:
$$
p_{k+1} - 2 p_k + p_{k-1} = (f^{tot}_k/m)  dt^2
$$

Assuming leap-frog integration scheme, we can substitute the velocities in the half step:
$$(v_{k+/2} - v_{k-/2}) dt = (f^{tot}_k/m)  dt^2$$
where we substituted:
 - $v_{k+/2} dt = p_{k+1} - p_k$
 - $v_{k-/2} dt = p_k - p_{k-1}$

To update velocity at half step by corrected position as $ v_{k+/2}  = (p_{k+1} - p_k) / dt$ seems to be valid approach. Nevertheless, we are not sure if the impulses introduced by constrain solver $dv_{k+/2}$ conserve linear and angular momentum.

$dv_{k+/2} = v'_{k+/2} - v_{k+/2} = (f^{int}_k/m)dt$
$dv_{k+/2} = (p_{k+1} - p_k)/dt -  v_{k-/2} - (f^{ext}_k/m)dt$
$dv_{k+/2} = (p_{k+1} - p_k + v_{k-/2}dt )/dt - (f^{ext}_k/m)dt$

## Projective dynamics as linear system

Projective dynamics minimize strain of constrains as well as action (Lagrangian) of the dynamical system. It is fomulated as linear system of equations $ A {\vec p} = {\vec b} $, comprising of stiffness of constrains $k_{ij}$ between points $p_i$ with mass $m_i$ and inertial term $I_i = m_i/dt^2$. The elements of matrix $A$ and vector $b$ can be written as follows: 

 - off-diagonal elements of PD marix: $A_{ij} = k_{ij}$
 - diagonal elements of PD marix: $ A_{ii} = I_i + \sum_j k_{ij}$
 - right hand side vector: $ b_i    = I_i p'_i + \sum_j{k_{ij} d_{ij}}$  
    , where:
     - $p'_i$ is the predicted position of particle $i$ by action of external forces $f^{ext}$ (i.e. ignoring the internal constraints $f^{con}$)
     - $d_{ij} = (p_j - p_i) (l_{ij}/|p_j - p_i|)$ is the optimal displacement of particle $i$ which would satisfy the constraint between particles $i$ and $j$ ( i.e. ignoring all other constraints).  

## Jacobi iteration as linear interpolation

The matrix equation $ A {\vec p} = {\vec b} $ of projective dynamics can be solved iteratively using Jacobi iteration scheme:

$$ 
p_i^{k+1} = \frac{ b_i -\sum_j A_{ij} p_j^{k} }{ A_{ii} } = \frac{ b_i -\sum_j k_{ij} p_j^{k} }{ I_i + \sum_j k_{ij} }
$$

If we update $b_i$ every iteration (which goes beyond the framework of solution of linear equations), then we can further expand this equation as follows: 

$$ 
p_i^{k+1} = \frac{ b_i^k -\sum_j k_{ij} p_j^{k} }{ I_i + \sum_j k_{ii} } 
          = \frac{ I_i p'_i + \sum_j{k_{ij} d^k_{ij}} -\sum_j k_{ij} p_j^k }{ I_i + \sum_j k_{ij} }
          = \frac{ I_i p'_i + \sum_j{k_{ij} ( d^k_{ij}} + p_j^k) }{ I_i + \sum_j k_{ij} }
          = \frac{ I_i p'_i + \sum_j{k_{ij} p'^k_{ij} } }{ I_i + \sum_j k_{ij} }
$$

Where we introduced $p'^k_{ij} = d^k_{ij} + p_j^k$ as the predicted position of particle $i$ which satisfy the constraint with particle $j$ ignoring all other constraints. Now if we ignore the inertia term $I_i \to 0$, it the above equation becomes simply weight average of predicted positions of particles due to all the constraints.
$$
p_i^{k+1} = \frac{  \sum_j{k_{ij} p'^k_{ij} } }{  \sum_j k_{ij} }
$$

To summarize, the Jacobi iteration of projective dynamics can be seen as weighted average of predicted positions of particles due to all the constraints weighted by stiffness and dynamically predicted positions weighted by inertia.

$$
p_i^{k+1} = \frac{ I_i p'_i + \sum_j{k_{ij} p'^k_{ij} } }{ I_i + \sum_j k_{ij} }
$$





