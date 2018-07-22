
# Object2Cell:

## Shared cell2object array:
 - **problem** - In typicall game-world we have:
    1. too many potential grid cells to store them all in single array
    2. multiple objects mapped into single grid cell
    3. single object is maped into multiple cells (sphere,boxes,lines,trinagle ... anything but point)
    All these issuse can be solved by use of `unordered_map` or implemented using either `chaing` (lists in buckets) or `open_adress_hashing` but it is not necessarily most efficient. Either we store arraylist in each bucket, or we have cells overlaping in array - both implementations have some problem https://stackoverflow.com/questions/31112852/how-stdunordered-map-is-implemented. If number of "find" calls dominates over "insert" calls we want better implementation. We would like construct hashmap<icell,{ioffset,nobjects}> pointing to array cellOrderedObjects[nobjects] so we can list all objects in array as follows:
    `for( i=0..cellmap[icell].nobjects ){ iobject = cellOrderedObjects+cellmap[icell].ioffset + i; }`
 - **Pseudocode point objects**
    - for all objects find icell and insert them to hashmap<icell,{ioffset,nobjects}>
        `for(oi : objects){ icell=tocell(oi); if(hashmap[icell]){ hashmap[icell].nobjects++; }else{ hashmap[icell].nobjects={ioffset=0,nobjects=1}; }; }`
    - iterate over hashmap and update ioffset
        `ioff = 0; for(cell : hashmap){ cell.ioffset=ioff; ioff+=cell.nobjects; cell.nobjects=0; }`
    - iterate over objects and store them to array cellOrderedObjects:```
            cellOrderedObjects = new int[nobjects];
            for(oi : objectts){ 
                icell=tocell(oi); 
                ioff =hashmap[icell].ioffset;
                cellOrderedObjects[ioff+hashmap[icell].nobjects] = oi
                hashmap[icell].nobjects+++;
            }
 - **Pseudocode object to multiple cells** - we simply generate fragments first
    - for all objects build fragments (part of object maped to single cells):
        `for(oi : objects){ 
            for( frag : tocells(oi) ){
                fragments.push({frag.icell,oi});
            }
        }`
    - for all fragments find icell and insert them to hashmap<icell,{ioffset,nobjects}>
        `for( frag : fragments){ icell=frag.icell; if(hashmap[icell]){ hashmap[icell].nobjects++; }else{ hashmap[icell].nobjects={ioffset=0,nobjects=1}; }; }`
    - iterate over hashmap and update ioffset
        `ioff = 0; for(cell : hashmap){ cell.ioffset=ioff; ioff+=cell.nobjects; cell.nobjects=0; }`
    - iterate over objects and store them to array cellOrderedObjects:```
            cellOrderedObjects = new int[nobjects];
            for( frag : fragments ){ 
                icell=frag.icell; 
                ioff =hashmap[icell].ioffset;
                cellOrderedObjects[ioff+hashmap[icell].nobjects] = frag.oi
                hashmap[icell].nobjects+++;
            }

# Projectile 2 Object
 - **Burst bounding box** - Burst of N projectiles is stored in shared bounding box (capsula? rather then bbox?)
 - **Time-step sweep** - Bounding box enclose both pos(t) and pos(t+dt) (i.e. position in subsequent time step)

## Direction sepatrated Sweep and Prune (SAP)
 - **rationale** - Since projectiles move fast in one direction, it's bounding box is prolonged along flight direction and quite narrow in perpendicular directions. Therefore Sweep and Prune (SAP) algorithm can be optimized if direction is selected perpendicular to flight direction. We therefore construct sevral disjuinct lists, each storing such bunches of projectieles for which projection of velocity in that direction is minimal.
 - **Pseudocode**
    - for all burstBox[j] find [i_min] with minimal direction[i].dot(burst[j].vel)
        - store burstBox[j] and its projection to direction[i_min] into burstList[i_min]
    - for all objectBox[j] calculate projection to all directions[i] and store to objectList[i]
    - sort all objectList[i] and burstList[i]
    - for all directions [i] sweep over burstList[i] and objectList[i] simultaneously and search collision
 - **choice of directions** - We would like to collide those projectiles with compact (sphere-like, box-like) objects which are scattered mostly in xz-plane ( vertical spread along y-axis is much smaller especially for objects on the ground). Also projectiles are often affected by gravity (especially bombs) and their trajectory is cureved in y-direction, thus v_y component is not conserved during flight. Threfore best choice of direction is azimuthal in xz-plane. We do not want to much directions as we have to keep sorted lists of all scene objects along each of them. Reasonable choice is therefore 4 directions ( x,z,zx,xy ) or 6-8 max.
- **Amortization of sorting cost** - due to temporal coherency (i.e. objects does not move so much between frames, especially for slow ground vehicles) the lists are almost sorted from previous frame. We use cheap updates by insertion sort, or shell sort (https://www.toptal.com/developers/sorting-algorithms/nearly-sorted-initial-order)



# AABB Tree
- https://www.azurefromthetrenches.com/introductory-guide-to-aabb-tree-collision-detection/
- https://github.com/lohedges/aabbcc

