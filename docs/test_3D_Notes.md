

## Game Sketches:

- test_BlockBuilder         - Can be a game about panel-building simulations - post apo game fighting in destructible socialistic panel building  

## Physical Problems

- test_CompressiveParticles - impact of compresible matrial on funel, usefull e.g. to simulate compression of nuclear material
- test_Electromagnetic      - plasma nozzle simulation, particles in electromagnetic field. TOKAMAK and other fuction reactor, Marnetic Nontainment 
- test_Radiosity            - full radiosity simulation in orthogonal coridor map, seems to work, but it is rather slow. worth looking how to optimize, then reuse for spacecraft.
- test_RayScattererMMC      - Monte-Carlo ray scattering in volumetric elements (tetrahedra), difinitely worth checking out
- test_Scatterer            - different scattering simulation, making all posible interaction between array of surface eleemnts (i.e. surface based, not volumetric)
- test_VortexLattice        - testing potential flow aerodynamics, around wing, can be very usefull to have some generator of aerocraft in javascript online on web. 

## Physical Problems (BAD)

- test_Elasticity           - simulation of girder like in bridge builder, but old way, direct force, inefficiet
- test_RigidBody            - rigid body simulation hanging on two strings, perhaps does not work correctly (inverted/trnasposred rotation matrix, or something)
- test_MultipoleAccel       - test Multiplole method, not sure if it work, not very interesting as a demo
- test_Stick                - dynamics of stick hitting object (kinda bsic)

## 3D geometry  (GOOD)

- test_AABBTree  - good demonstration of effificne spatial indexing
- test_Mesh      - various mesh operation, not sure which, probably old, it would be good to reinvestigate
- test_Patches   - smooth triangular patches, look how this is implementd, I fogrot
- test_SphereSampling - mapping 2D map of sphere using ocahedro and icosahedra, very useful for globus-like maping, plants, asteroides, but also ray-scattering problem (maping histrogram of directional incident ray intensity)
- test_SphereTree  - Diffusion Limited Aggregation in 3D, grow interesting fractal structure which are accessible from outside, yet connected. Alternative to Hydraulic Erosion for terain generation. Can be easily made miulti respolution (low res DLA, high-res hydraulic erosion) 

## 3D geometry (BAD)

- test_MusculeEditor     - Editor of muscule like shapes / like Loft between curves. Not sure how usefull. Maybe good for organic shapes,  
- test_SphereGaussSeidel - packing boxes 

## Less Interesting / Techinal BAD

- test_Camera          - camera reycasting and picking
- test_Collision       - ridig body on terrain - it does not work properly. Perhaps rigid body simulation is somehow inverted (transposed rotation matrix or somethig?)
- test_EditorGizmo     - How to use manipulation Gizmo
- test_GUI             - How to use GUI
- test_MousePicking    - seems just like test_Camera, picking particles by mouse
- test_Projection      - this is just for debugging the camera view frusturm
- test_QuatRotSampling - debugging how to aling quaternion, probably some old method, now we know more efficient ways
- test_Raytracing      - picking trinagles by a ray from camera 
- test_Solids          - Solids.h, CMesh.h demonstration
- test_SceneGraph      - SceneGraph for modeling (duplicate gemoetry etc.)
- test_TrussBuilder    - truss builder similar to test_BlockBuilder




