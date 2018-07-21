
## non-Rotational BBoxes 

### Spherical BBox
 - minimal storage requirement (just 4 numbers), fits in vec4/float4 on GPU
    - Since most objects already has `pos`, it needs to add only one additional scalrar `r`
 - projection to all diractions is the same, therefore we do not have to care about rotation of object neigher sphere
 - Sphere-sweep is simple and efficient ( line-point resp. line-line distance )
 - prolonged objects ( e.g. gun,rocket, ship, train vagoon ) cannot be tightly enclose leading to efficiency issues

### Axis Alligned Bounding Box (AABB)
 - good for storing large static chunks of scene
 - needs to be updated when object rotated
 - evaluating projections along non-axis-alligned is relatively costly ( 8-points )

## Rotation-dependent BBoxes

### Cylinder & Capsula
 - rotational structure with minimal storage requrements (just one direction vector instead of Matrix)
    - since objects already has position & direction, only two additional numbers (radius,lenght) are required
 - can stroe effienctly both proponged objects (shot burst) and flatened objects ( disk, city, panel )
 - Capsula vs. Cylinder
    - capsula enclosing n-points is typically more compact
    - capsula-sphere overlap is eariser to calculate as line-point distance
    - capsula is more difficult to construct than cylinder ( needs sqrt per each point     `tc = sqrt(Rc^2-ri^2) - ti`      )

### Rotated BBox
 - Often most thightly-fitting shape from all easy-to-calculate shapes
 - Can rotate easily with objects
 - evaluating overlap of two differently rotated shapes may be rather complicated and costly

### Elipsoiede
 - it is just sphere in distorted space
 - however interaction of two differently rotated elipsoides can be difficult - it can be transformed to sphere-elipsoide intersection
