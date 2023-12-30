# Particle In Cell simulation

Particle in cell is class of algorithms which approximates particle-to-particle interactions by particle-to-grid interactions. This is efficient because while particle-particle interactions scales $O(n^2)$ with number of particles $n$, each particle interacts just with constant number of grid-points therefore reducinc complexity to $O(n)$. Typically the grid is rectangular, which means that cell are boxes with $2^D$ corners (e.g. 4 in 2D, 8 in 3D).


## Algorithms

* **Projection** - we map the particles into corresponding boxes. We update particle densities, mass densities, charge densities and momentum densities
* **Grid-Field solver**
* **Interpolation**


## Collision of particle to grid



### Momentum Conservation






* [Particle in Cell in shader](https://www.shadertoy.com/view/XcB3zm)
* [Point Distance Function](https://www.shadertoy.com/view/XcSGzm)