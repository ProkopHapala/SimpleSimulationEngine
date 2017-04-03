
# TO DO

 - global search for 1 molecule on rigid substrate
    -  on rigid substrate we can efficiently use fast methods for collision detection like
        - GridFF distance filed approximation
        - GridMap accelerated short-range NBody  
    - Do this using OpenCL (fully inside GPU) that involves:
        - GridFF evaluation, GridFF sampling, RigidBody dynamics, 
        - First is evaluated just pauli repulsion (Coarse grid), than 
        -   

## Fast Forcefield

#### Shoothened fields
replacing $1/r$ by Lorenzian $1/\sqrt{r^2 + w_0^2}$ not only removes singularity but it is actually more realistic for atoms which are not point-charges but clouds of electron density.

#### Shoothened fields
$(1-1/r^2)^2-1$ - looks very similar to LJ but can by simply evaluated from $1/r^2-1$ sored in GridFF. 


### GridFF for rigid substrate

This will reduce NBody problem with $O(N^2)$ complexity to simple grid sampling. However grid based forcefiled has also several problems. List of pro-vs-cons:

 - Pros:
    - force for each particle can be evaluated in constant time, even for long-range forces (e.g. electrostatics)
 - cons:
    -  cannot be very fine since memory it scales as $M(n_x.n_y.n_z)$ with number of grid samples
    -  non-linear properties of atoms can be hardly encoded
        - electrostatic potneital $E_{el}(r)=q_iq_j/r$ can be easily encoded. We simpley store $V(r)=\sum_i q_i/|r-r_i|$ 
           
        - but Lenard-Jones potential $E_{LJ}(r)=e_{ij}((R_{ij}/r)^{12}-2(R_{ij}/r)^6)$ cannot be easily encoded 

**distance field approximation** (DFA-FF) for collision detection exploit the fact that short range potentials such as **Pauli Repulsion** decay fast with distance. Therefore good approximation of total force is given jus by considering neares atom. We usually just want to encode the boundary of vdW radius of atoms $R_i$. We should store $V(r) = \sum_i 1/(|r-r_i|^2 - R_i^2)$ or $V(r) = \min_i 1/(|r-r_i|^2 - R_i^2)$. Since this is very crude approximation it is usefull just for fast pre-screening we should use very coarse grid. Grid step should be just ~0.25-0.5 A (which leads to 8x-64x lower memory consumption than typical step 0.1A).

### Few-large chunks

GridFF and DFA-FF can be used not only for acceleration of rigid substrate but also for interaction of two large rigid molecules (e.g. two nano-crystals, viruses etc.). Atoms of one interacts with GridFF of the other and vice-versa. It is also possible to calculate overlap of the two grids directly, but that is probably slower, because there is much more grid-points than atoms. Also precission of grid-grid interaction would be lower than for grid-atom interaction. Since grids are approximations of atom-wise interactions.
