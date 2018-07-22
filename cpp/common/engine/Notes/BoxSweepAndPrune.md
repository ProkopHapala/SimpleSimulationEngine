
## Motivation

Sweep-and-prune algorithm sorts objects just along one axis (say x-axis), therefore there can be still many objects seemingly overlaping (since they have similar x-component), even though they are very far along other axis (y,z). This is bad because:
    - all m objects overlaping on sweep-axis pairwise ( O(m^2) ) if they overlap also along y,z;
    - whole sweep-axis has to be sorted at-once ( O(NlogN))

## Box-Sweep-Prune with const-length boxes
 - We can first put objects into smaller 1D cells of fixed size along y-axis  (O(N) cost), and then sort content of each box just along x-axis
 - If we wan to use fixed size array to store permutations of particles we may do it in similar way as in case ov 3D-grid method

## Box-Sweep-Prune with const-number boxes
 - We may also make cells such that each contains fixed number of particles., therefore we can easily store permutation in fixed-size array
    - however this means we have to sort all objects along y-axis first, therefore we loose advantage for ommition of sorting
