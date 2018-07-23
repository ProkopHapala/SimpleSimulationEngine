
## Fully Branching AABB Tree (each leaf has 2 forks)

## Leaf branch Distiction
 - top bit of index can denote if branch is Leaf or Node
    - this can be also +/- sign if we use signed int

## Array Storage
 - We can gurantie heap condition ( for each node A, all its sub-nodes have higher id than A ) therefore we can store sub nodes in dynamic array
 - We can also distribute members (bbox, index, pointer) to separate arrays
 - We can also add additional proeprties which are rarely used ( e.g. volume/area/cost of given bbox ) into axuliary array

## cache consideration - Binarry vs N-tree
 - If branches are stored as pointers (resp.indexes) than we do not gain much better cacheing with N-treethan with binary tree, since we have to dereference all N bboxes anyway.
 - Much better cacheing can be achieved if we point to beginning of all N bboxes
    - If we use fully-branchind nodes (i.e. each node is guarantied to store maximum number (N) branches) than we are fine
    - If we can in principle store number of bboxes in block, but that will cause problems with allocation of new boxes

## Balancing
 - It is rather easy to balance by swaping sub-nodes on the same level, especially for binary tree
    Consider 2 levels with 4 nodes, there are 3 possible arrangements:
        AA  AB  AB
        BB  AB  BA
    - therefore given configuration ((1,2),(3,4)) we generate new ((1,3),(2,4)) and ((1,4),(2,3)) and test if cost of such configuration is lower.


