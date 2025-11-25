## User 1

We need to make some high-level code review of our SpaceCraft building system. 

Our goal is to make comprehensive high-level overview of current state of the system, and what needs to be done, and where we find relevant functionality. This should be all sumarized in new markdown document SpaceCrafting_new.md

Our primary focus now is on bulding the spacecraft - that if creating the Mesh from Lua script 
using @EditSpaceCraft.h , @SpaceCraft.h , @SpaceCraftComponents.h 
which is later exported into Mesh using @MeshBuilder2.h and only then it may be exported to truss for simulation

 (which is then exported to Truss for simulation but that is less relevant), also we will later simulated both dynamics using 
@TrussDynamics_d.h and Ridiative Heat transfer and ratiation scattering using @TriangleRayTracer.h and @Radiosity.h , but these are less relevant taks for now (we should only note this as background, context)

There are some documents already made which you should read and note in SpaceCrafting_new.md
you should extract most relevant information and sumarize it there (before examining actuall code)

The usefull .md files are:

@spaceCraft_mesh_export_cli.md 
@SpaceCrafting.md @SpaceCraftingWithBlocks_new.md @TriangleOcclusionRaytracer.md 
@ConstructionBlock.md @MeshBuilder2.md @ScatteringSolverComparison.md 

Now the code which use and run it (the build targets, entry points) are 
@constructionBlockApp.cpp @spaceCraftDynamics.cpp @spaceCraftEditor.cpp @SpaceCraftEditorNew.cpp @spaceCraftMeshExport.cpp 

please make intial overview and then I will give you more specific what to describe in the documentation

---

## User 2

This is perfect ! Thank you ! You nailed it !
You even captured many things on which I wanted you to focus. 

As you pointed out there are several open chalanges on which we should focus now. In what state of implementation and resolution these are? what are the tools in our disposal to solve them? We should investigate this deeper.

A) Lua mapping kinda works, we do not need to go into more details.

B) yes this I wanted to stress and investigate further what is the relation and differences between the older spaceship generator (Build_truss or so) and the new Bock-based generator.
As far as I know I originally had point-like nodes (single mass-point) and the girder were pointy (they started and ended in that single vertex), the motivation was to make it simple to edit (you can take any truss (i.e. points connected by lines) and generate more detailed truss by replacing the lines by girders. However this had a serious problem: The so generates spacecraft was very "wobbly" because the single-point nodes does not constrain or tranfer any bending or trosion force (constrain only distance). I was thinking that rope will be enough, but were not, the mechanis of the ship was unstable, difficult to converge in elesticity simulator, dificult to aim, dificult to damp oscillation. Therefore I decided to rework the whole system using solids (like cube, or octahedron) as the nodes for the truss. Now the girder are attached to multiple points, for example to quad face of cube, which gives it more stability, especially in bending and twisting (torision), efficiently damping the vibration modes. It is much more stable now. but I did not finished yet the problem how to properly attach all other compnent to the new system. 
1) where exactly attach ropes?
2) where exactly attach slides?
3) where exatly attach plates (radiators, shields)

C) yes, exactly, the sliders / paths are the core issue, we always have problem with these and we need to carefully document how they are done, and also how they are connected to the actuall dynamical simulation in @TrussDynamics_d.h , notice we have these axuliary EdgeVertBond
TrussDynamics_d::evalEdgeVerts(){ which we aplly in 
TrussDynamics_d::run_LinSolve(int niter) {
as external force, which are therefore not part of the Truss, but are considered as some external soft bonds constrains which can move. Proper implementation of them is crucial. they move alnong the paths - where the paths are currently just selection of points on the truss, and are linearly interpolated (therefore these are actually edges between the truss vertexes) - a polyline. 
... perhas to make the movement more smoother it would be good to replace the linar interpolation by e.g. cubic splines (Hermite or Bspline). When we export .truss file in @spaceCraftMeshExport.cpp we should note we need to export also these axuliary links / drivers, otherwise the truss is incomplete. For example typically wheels are attached to girders only using these soft connections. So if we ommit them the whole spacecraft will fall appart (wheels will dedatch). We also want to have general system which allows us to attach not only wheels but also e.g. ropes or other girdes to the paths defined on girders and wheels.
One of the chalanges we had, was that the there is some system to "find side" on the girder, which select the closest edge-strip when building the script forom Lua of from editor. We need to look on this in detail while debugging because I remember it was notoriously buggy, and we spend lot of tim on it, and I'm still not sure it works.

D) we do not need to care about this. It is not importaint

 I want to add two more points:

E) Plates (Shields, Radiators) - I'm not sure in which phase of implementation we are, we perhaps only started. There are two questions/chalanges. 
1) How to attach these to girders (or eventually to ropes), general idea of the high-level script was to set some points along the lengh of girder as starting and ending of the plates. And then to generat triangular grid with some sampling to fill it up. However this either means we need to match the subdivision of the girdrs (having same number of segments on each side) or how else we can ensure the teselation of the plat match segments of the girder on both sides ?
2) should we offset (inset) the girders and then connect them by some axuliary line segments (sticks). Maybe we can implement both options using two different mesh-generating functions. 
3) should we give these plates some structura stiffnes in bending, which means non-zero thickness? For that purpose we have created some slab-generators which we tested in 
@constructionBlockApp.cpp and in @DrawUV.h , we should note it. For some components adding flexular stifnes may be overkill (like for radiuators or solar-sails, but e.g. for shields, or pusher plates, it may be essential.

F) Nozzles / Dishes

I'm not sure if we already started proper implementation of trusters which are typically made of quite large pusher-plate (for Orion-prokekt like nuclear spacecraft (which is basically sturdy ridig shield), or by magnetic nozzle for project Daedalus type of spacecraft with magneto-interial containment of nucleary pulse plasma). Both pusher plat and  magnetic nozzle should have shape of disk or parabola and can be generated again using  function tested in @constructionBlockApp.cpp and implemented in @DrawUV.h . They should be generally more rigid than the shields or radiators, so they should definitely have some thickness (double-layer slab conected by cross lineks)

Perhas the biggest chalange with these structures is how to connect them to the rest of the ship. 

On one hand to make the spacecraft solid (not wobly) we should connect it firlmy at mitliple points so it cannot bend or twist. 

On the other hand, ofthen we want to have some dynamical damping ( to damp the impulses from nuclear pulse engine), these pulses have very different range of amplitudes, and therefore different damping lengh  ( hundrest of meters for orion-like fission bonds, and only few menters or centimeters for micro-fussion pelets, some magnetic nozzles may be fully continuous and require no damping. ). For large amplitude damping we should attach the ship by ropes. For small amplitude damping we can direcly attach them byt solid sticks (normal rigid edges of the truss) but made from special damping material (e.g. instead of kevlar, or steel we select "dampler" for the stic in the     SpaceCraftWorkshop

But it is also about the construction - should we make these connections manually one-by one (connecting individual vertexes of the thruster nozzle with individual vertexes of the girder ? or should we use some automatic generator. I thin in the @MeshBuilder2.h we have some "weld" function which inside some box find all vertexes between the two selections of vertices (from one and other compnent) and connect all with all (under certain distance trashold). If not we should implemnt this. Anyway this is rather crude approach. But manula connection is laborious. 



