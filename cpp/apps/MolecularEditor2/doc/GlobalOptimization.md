
## Macro State Tree

 - Configuration space is splited into disjunct region acording of subsed of descriptors. 
 - Bins are stored in some spatial indexing datastructure (typically hasmap) which does not require allocate non-sampled bins.
 - For each sampled tile object is allocated which contains
    - Free energy of the bin (the higher variability and the lower energy, the more attractive is the region)
    - Probability to hit configuration which was already sampled (this lowers the attractivity)
 - During sampling regions are selected by some algoritm, like Monte Carlo according these two marcoparameters



