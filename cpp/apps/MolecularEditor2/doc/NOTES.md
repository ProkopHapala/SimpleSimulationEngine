
### TODO:

 * Try conf dynamics with single $H_2O$ molecule on rigid NaCl substrate
 * checks
    * check angular forces (not acos)
    * check size of forces

#### Rigid Body:

Rigid bodies are integrated with soft-body molecules in following way:
 * There is common array of force and position `apos` and `aforce`
 * coordinates and forces are transfromed by functions:
   * `frags2atoms()` - generate atomic coordinates `apos` from molecular `poses`
   * `aforce2frags()` - project atomic forces `aforce` to `poseFs`
 * Thanks to this atomic froces `aforce` can be computed by common function e.g. `eval_MorseQ_On2` and `eval_FFgrid` for both *rigid* and *flexible* molecules
 * Only problem is that optimizer `DynamicOpt` now have to optimize both *atomic* and *molecular* coordinates. Some parts of `apos` should be excluded from `DynamicOpt`. 
    * If the arrays are organized like  `|..rigid..|..flexible..|` it can be solved by simple shift of the pointer or modified lenght. Even better if we use common buffer organized like `|..pose..|..rigid..|..flexible..|` 


#### Grid potential:

There is ractorized (separable) Morse potential on a grid. Common array of force is interpolated and only than multiplied by coefficients correspoding to particular atom type

Grid can easily accelerate any additive potential from atoms of rigid substrate and single particle. But if there are different types of particles than acceleration of Lennard-Jones potential may be problmatic sice it's dependence on parameters of atoms (in particular radius $R_i$ ) is inhomogenous. 

[Coulomb](https://en.wikipedia.org/wiki/Coulomb%27s_law)
$$E_a(r) = \sum_i q_a q_i /r$$

[Lenard-Jones](https://en.wikipedia.org/wiki/Lennard-Jones_potential)
$$E_a(r) = \sum_i e_a e_i (  (R_a + R_i)^{12}/r^{12} - 2(R_a + R_i)^6/r^6 )$$

in this case $R_a$ cannot be easily factorized in front of the expression
[Mores](https://en.wikipedia.org/wiki/Morse_potential)

$$E_a(r) = \sum_i e_a e_i ( \exp(-2a(r-R_a-R_i) - 2\exp(-a(r-R_a-R_i) )$$
factorization:

$$E_a(r) = \sum_i e_a e_i  ( \exp(-2a(r-R_i)\exp(2aR_a) - 2\exp(-a(r-R_i)\exp(aR_a) )$$ so it can be easily expressed as sum of two factorized potentials:

$$E^+_a(r) = e_a \exp(2aR_a) \sum_i  e_i \exp(-2a(r-R_i) ) $$
$$E^-_a(r) =  e_a exp(aR_a)  \sum_i 2e_i  \exp(-a(r-R_i) )$$


