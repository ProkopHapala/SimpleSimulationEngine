
#ifndef  SpaceCraft2Mesh2_h
#define  SpaceCraft2Mesh2_h

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "Noise.h"
//#include "SphereSampling.h"
//#include "DrawSphereMap.h"

//#include "Truss.h"
#include "MeshBuilder2.h"
//#include "MeshBuilderDraw.h"
#include "SpaceCraft.h"
#include "TrussDynamics_f.h"
#include "TrussDynamics_d.h"
#include "testUtils.h"

#include "Solids.h"
using namespace SpaceCrafting;

namespace Mesh{

static double g_Rcollapse =   0.1;
static double g_Ranchor   =  10.0;
static int    g_AnchorStickType =  0;

void exportSim( TrussDynamics_d& sim, const Builder2& mesh, const SpaceCraftWorkshop& shop, int nneighmax_min=0 ){
    // {
    //     printf("#====== Materials \n");      for(const Material*      o: shop.materials.vec     ){  o->print();    }
    //     //printf("#====== PanelMaterials \n"); for(const PanelMaterial* o: shop.panelMaterials.vec){  o->print();    }
    //     printf("#====== StickMaterials \n"); for(const StickMaterial* o: shop.stickMaterials.vec){  o->print();  }
    //     //exit(0);
    // }
    int np = mesh.verts.size();
    int nb = mesh.edges.size();
    printf( "exportSim() START nvert=%i nedge=%i \n", np, nb );
    mesh.checkAllPointsConnected(true, true);
    int* nneighs = new int[ np ];
    //printf( "exportSim() np=%i\n", np );
    // find max number of neighbors
    for(int i=0; i<np; i++){ nneighs[i]=0; }
    for(int i=0; i<nb; i++){ const Vec2i& e = mesh.edges[i].lo; nneighs[e.a]++; nneighs[e.b]++; }

    int nneighmax = nneighmax_min;
    for(int i=0; i<np; i++){ int ni=nneighs[i]; if(ni>nneighmax)nneighmax=ni; }

    int nneigh_hist[nneighmax+1];
    for(int i=0; i<=nneighmax; i++){ nneigh_hist[i]=0; }
    for(int i=0; i<np; i++){ nneigh_hist[ nneighs[i] ]++; }
    //for(int i=0; i<=nneighmax; i++){ printf( "exportSim() nneigh_hist[%i] %i \n", i, nneigh_hist[i] ); }
    //exit(0);

    //printf( "exportSim() nneighmax %i \n", nneighmax );
    sim.recalloc(np, nneighmax, nb );

    // --- EdgeVerts
    //printf( "exportSim() sim.nEdgeVert %i \n", sim.nEdgeVert );


    Vec2d kPullRange  = Vec2d{1.e+300,-1.e+300};
    Vec2d kPushRange  = Vec2d{1.e+300,-1.e+300};
    Vec2d l0Range     = Vec2d{1.e+300,-1.e+300};
    Vec2d massRange   = Vec2d{1.e+300,-1.e+300};
    std::unordered_set<int> stickTypes;;

    // fill neighs
    for(int i=0; i<np; i++){ sim.points[i].f=mesh.verts[i].pos; sim.points[i].e=0.0f; }
    for(int i=0; i<np; i++){ nneighs[i]=0;   }
    for(int i=0; i<sim.nNeighTot; i++){ sim.neighs[i]=-1; }
    if(sim.neighBs ){ for(int i=0; i<sim.nNeighTot; i++){ sim.neighBs [i]=(int2){-1,-1}; } }
    if(sim.neighB2s){ for(int i=0; i<sim.nNeighTot; i++){ sim.neighB2s[i]=0;             } }
    for(int i=0; i<nb; i++){
        //printf( "exportSim()[%i] \n", i );
        const Quat4i& e = mesh.edges[i];
        int ia = e.x*sim.nNeighMax + nneighs[e.x];
        int ib = e.y*sim.nNeighMax + nneighs[e.y];
        sim.neighs[ ia ] = e.y; 
        sim.neighs[ ib ] = e.x;
        if(sim.neighBs){
            sim.neighBs[ ia ] = (int2){e.y,i}; 
            sim.neighBs[ ib ] = (int2){e.x,i};
        }
        if(sim.neighB2s){
            sim.neighB2s[ ia ] =  (i+1); 
            sim.neighB2s[ ib ] = -(i+1);
        }

        //printf( "e.w %i \n", e.w );
        if(e.w<0){ printf( "ERROR in exportSim() mesh.edges[%i].type=%i \n", i, e.w ); exit(0); }
        if(e.w>=shop.stickMaterials.vec.size()){ printf( "ERROR in exportSim() mesh.edges[%i].type=%i > stickMaterials.size(%i)\n", i, e.w, e.w>=shop.stickMaterials.vec.size() ); exit(0); }
        const StickMaterial& mat = *shop.stickMaterials.vec[e.w];
        stickTypes.insert(e.w);
        // l0, kPress, kPull, damping
        //double l0 = (mesh.verts[e.y].pos - mesh.verts[e.x].pos ).norm();
        double l0 = (sim.points[e.y].f - sim.points[e.x].f ).norm() * ( 1.f - mat.preStrain );
        double mass = l0*mat.linearDensity;
        sim.points[e.x].w += mass*0.5;
        sim.points[e.y].w += mass*0.5;
        Quat4d param = (Quat4d){ l0, mat.Kpush/l0, mat.Kpull/l0, mat.damping }; 

        if(isnan(param)){ printf("ERROR in exportSim() mesh.edges[%i].type=%i param(l0=%g, kPush=%g, kPull=%g, damping=%g) \n", i, e.w, param.x, param.y, param.z, param.w ); exit(0); }

        // if mat.name begin with "rope"  
        if( strncmp( mat.name, "rope", 4 )==0 ){
            sim.damped_bonds.push_back( i );
            sim.damped_points.insert( e.x );
            sim.damped_points.insert( e.y );
        }
        l0Range   .enclose( param.x );
        kPushRange.enclose( param.y );
        kPullRange.enclose( param.z );
        massRange .enclose( mass    );

        //param.z = 1e+7;
        //if( param.z>1e+7) param.z=1e+7;   // This Makes it converge
        if( (param.z<1e+6) || (param.z>1e+8) )
        {
            //printf( "exportSim ib: %6i typ: %2i %-10s l0[m]: %6.2f m[kg]: %6.2f k[N/m]: %.3e %.3e \n", i, e.w, mat.name, mass, param.x, param.y, param.z, param.w );
            //printf( "exportSim ib: %6i typ: %2i %-10s   length: %6.2f [m] mass: %6.2f [kg] kPush: %.3e [kN/\%] kPull: %.3e [kN/\%] \n", i,  e.w, mat.name,  l0, mass, mat.Kpush*1e-6, mat.Kpull*1e-6 );
            //( "exportSim ib: %6i typ: %2i %-10s   length: %6.2f [m] mass: %6.2f [kg] m/l %6.2f [kg/m] k: %.3e [N/m] k/l: %.3e [10kN/\%] \n", i,  e.w, mat.name, l0, mass, mat.linearDensity, param.z, mat.Kpull*1e-6 );
        //    const Material& M = *shop.materials.vec[mat.materialId];
        //    printf( "stick[%i] par(%7.3f,%5.2e,%5.2e,%5.2e) Stick(%s,%g[m^2],%g[m])K(%5.2e,%5.2e) Stick()mat(%s,K(%5.2e,%5.2e))\n",  i, param.x, param.y, param.z, param.w,   mat.name, mat.area, mat.diameter, mat.Kpush, mat.Kpull, M.name, M.Kpull, M.Kpush );
        }
        sim.params[ ia ] = param;
        sim.params[ ib ] = param;
        nneighs[e.x]++; nneighs[e.y]++;
        if(sim.bonds){
            sim.bonds[i]     = *(int2*)&e.lo;
            sim.bparams[i]   = param;
            //sim.l0s[i]       =  l0;
            //if(i==6272){ printf( "exportSim [ib=%i](%i,%i) param.x=%g l0=%g \n", i, e.x,e.y, param.x, l0 ); }
            sim.maxStrain[i] = (Vec2d){ (mat.Spull/mat.Kpull), (mat.Spush/mat.Kpush) };
            sim.strain[i] = 0;
        }
    }

    sim.checkMasses();

    sim.cleanForce();
    sim.cleanVel();
    //for(int i=0; i<sim.nPoint; i++){ sim.points[i].f.addRandomCube(0.1); }

    printf( "exportSim() StickMaterials: \n" );
    for(int i: stickTypes){ shop.stickMaterials.vec[i]->print(); }
    printf( "exportSim() kPushRange %g %g \n", kPushRange.x, kPushRange.y );
    printf( "exportSim() kPullRange %g %g \n", kPullRange.x, kPullRange.y );
    printf( "exportSim() l0Range    %g %g \n", l0Range.x,    l0Range.y    );
    printf( "exportSim() massRange  %g %g \n", massRange.x,  massRange.y  );
    printf( "exportSim() DONE! \n" );
    delete [] nneighs;
    //exit(0);
}


void exportSim( TrussDynamics_f& sim, const Builder2& mesh, const SpaceCraftWorkshop& shop ){
    // {
    //     printf("#====== Materials \n");      for(const Material*      o: shop.materials.vec     ){  o->print();    }
    //     //printf("#====== PanelMaterials \n"); for(const PanelMaterial* o: shop.panelMaterials.vec){  o->print();    }
    //     printf("#====== StickMaterials \n"); for(const StickMaterial* o: shop.stickMaterials.vec){  o->print();  }
    //     //exit(0);
    // }
    int np = mesh.verts.size();
    int nb = mesh.edges.size();
    printf( "exportSim() START nvert=%i nedge=%i \n", np, nb );
    int* nneighs = new int[ np ];
    //printf( "exportSim() np=%i\n", np );
    // find max number of neighbors
    for(int i=0; i<np; i++){ nneighs[i]=0; }
    for(int i=0; i<nb; i++){ const Vec2i& e = mesh.edges[i].lo; nneighs[e.a]++; nneighs[e.b]++; }
    int nneighmax = 0;
    for(int i=0; i<np; i++){ int ni=nneighs[i]; if(ni>nneighmax)nneighmax=ni; }

    int nneigh_hist[nneighmax+1];
    for(int i=0; i<=nneighmax; i++){ nneigh_hist[i]=0; }
    for(int i=0; i<np; i++){ nneigh_hist[ nneighs[i] ]++; }
    //for(int i=0; i<=nneighmax; i++){ printf( "exportSim() nneigh_hist[%i] %i \n", i, nneigh_hist[i] ); }
    //exit(0);

    //printf( "exportSim() nneighmax %i \n", nneighmax );
    sim.recalloc(np, nneighmax, nb );

    // --- EdgeVerts
    //printf( "exportSim() sim.nEdgeVert %i \n", sim.nEdgeVert );

    // fill neighs
    for(int i=0; i<np; i++){ sim.points[i].f=(Vec3f)mesh.verts[i].pos; sim.points[i].e=0.0f; }
    for(int i=0; i<np; i++){ nneighs[i]=0;   }
    for(int i=0; i<sim.nNeighTot; i++){ sim.neighs[i]=-1; }
    if(sim.neighBs ){ for(int i=0; i<sim.nNeighTot; i++){ sim.neighBs [i]=(int2){-1,-1}; } }
    if(sim.neighB2s){ for(int i=0; i<sim.nNeighTot; i++){ sim.neighB2s[i]=0;             } }
    for(int i=0; i<nb; i++){
        //printf( "exportSim()[%i] \n", i );
        const Quat4i& e = mesh.edges[i];
        int ia = e.x*sim.nNeighMax + nneighs[e.x];
        int ib = e.y*sim.nNeighMax + nneighs[e.y];
        sim.neighs[ ia ] = e.y; 
        sim.neighs[ ib ] = e.x;
        if(sim.neighBs){
            sim.neighBs[ ia ] = (int2){e.y,i}; 
            sim.neighBs[ ib ] = (int2){e.x,i};
        }
        if(sim.neighB2s){
            sim.neighB2s[ ia ] =  (i+1); 
            sim.neighB2s[ ib ] = -(i+1);
        }
        
        //printf( "e.w %i \n", e.w );
        if(e.w<0){ printf( "ERROR in exportSim() mesh.edges[%i].type=%i \n", i, e.w ); exit(0); }
        if(e.w>=shop.stickMaterials.vec.size()){ printf( "ERROR in exportSim() mesh.edges[%i].type=%i > stickMaterials.size()\n", i, e.w, e.w>=shop.stickMaterials.vec.size() ); exit(0); }
        const StickMaterial& mat = *shop.stickMaterials.vec[e.w];


        // l0, kPress, kPull, damping
        //double l0 = (mesh.verts[e.y].pos - mesh.verts[e.x].pos ).norm();
        double l0 = (sim.points[e.y].f - sim.points[e.x].f ).norm() * ( 1.f - mat.preStrain );
        double mass = l0*mat.linearDensity;
        sim.points[e.x].w += mass*0.5;
        sim.points[e.y].w += mass*0.5; 
        Quat4f param = (Quat4f){ l0, mat.Kpush/l0, mat.Kpull/l0, mat.damping }; 
        {
        //    printf( "exportSim(ib=%i) length=%g[m] mass=%g[kg] fPush(%g[10kN~ton]\%1) fPull(%g[10kN~ton]\%1)\n", i, l0, mass, mat.Kpush*1e-6, mat.Kpull*1e-6 );
        //    const Material& M = *shop.materials.vec[mat.materialId];
        //    printf( "stick[%i] par(%7.3f,%5.2e,%5.2e,%5.2e) Stick(%s,%g[m^2],%g[m])K(%5.2e,%5.2e) Stick()mat(%s,K(%5.2e,%5.2e))\n",  i, param.x, param.y, param.z, param.w,   mat.name, mat.area, mat.diameter, mat.Kpush, mat.Kpull, M.name, M.Kpull, M.Kpush );
        }
        sim.params[ ia ] = param;
        sim.params[ ib ] = param;
        nneighs[e.x]++; nneighs[e.y]++;
        if(sim.bonds){
            sim.bonds[i]     = *(int2*)&e.lo;
            sim.bparams[i]   = param;
            //sim.l0s[i]       =  l0;
            //if(i==6272){ printf( "exportSim [ib=%i](%i,%i) param.x=%g l0=%g \n", i, e.x,e.y, param.x, l0 ); }
            sim.maxStrain[i] = (Vec2f){ (float)(mat.Spull/mat.Kpull), float(mat.Spush/mat.Kpush) };
            sim.strain[i] = 0;
        }
    }
    sim.cleanForce();
    sim.cleanVel();
    //for(int i=0; i<sim.nPoint; i++){ sim.points[i].f.addRandomCube(0.1); }
    printf( "exportSim_f()  DONE! \n" );
    delete [] nneighs;
    //exit(0);
}

int exportBuckets( const SpaceCraft& craft, Buckets* buckets=0, int nPerBucket=16, bool bHardLimit=true ){
    // Goes over all the StructuralComponents in the SpaceCraft and exports them to the Buckets (i.e. bounding boxes), split each structural component into a number of buckets.
    // if no output Buckets is specified, it just counts the number of buckets
    //printf("exportBuckets(out=%i,nPerBucket=%i,hardlimit=%i) nGirder=%i nRings=%i nRopes=%i \n", buckets!=0, nPerBucket, bHardLimit, craft.girders.size(), craft.rings.size(), craft.ropes.size() );
    int nBuck = 1; // bucket 0 is reserved for un-assigned objects
    for(Girder* o: craft.girders){ nBuck += o->toBuckets( nBuck, nPerBucket, buckets, bHardLimit ); }
    for(Ring*   o: craft.rings  ){ nBuck += o->toBuckets( nBuck, nPerBucket, buckets, bHardLimit ); }
    for(Rope*   o: craft.ropes  ){ nBuck += o->toBuckets( nBuck, nPerBucket, buckets, bHardLimit ); }
    return nBuck;
}

void applySliders2sim( SpaceCraft& craft, TrussDynamics_d& sim, double* control_speed ){
    double dt = sim.dt;
    for( int i=0; i<craft.sliders.size(); i++ ){
        Slider* o = craft.sliders[i];
        int icon = o->icontrol; if(icon>=0){  o->speed = control_speed[icon]; }  // Controls
        EdgeVertBond& ev = sim.edgeVertBonds[i];
        Vec3d  d  = sim.points[ev.verts.y].f - sim.points[ev.verts.x].f;
        Vec3d  dv = sim.vel   [ev.verts.z].f - sim.vel   [ev.verts.y].f*ev.c + sim.vel[ev.verts.x].f*(1-ev.c);
        double l  = d.norm();
        double f  = d.dot( ev.f )/l; // force along the slider path
        double v  = d.dot( dv   )/l; // velocity along the slider path
        o->move( dt, l, v, f );
        // -- update corresponding EdgeVerts
        ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );
        ev.verts.z = o->ivert;
        //ev.K = 1000.0;
        //o->updateEdgeVerts( sim.points );
    }
}

void sliders2edgeverts( SpaceCraft& craft, TrussDynamics_d& sim ){
    // NOTE: Slider paths must be initialized before calling this function.
    // This is typically done by calling `o->updateSlidersPaths(...)` for each Ring.
    // ----- Sliders to EdgeVerts
    sim.nEdgeVertBonds = craft.sliders.size();
    if(sim.nEdgeVertBonds == 0) return;
    sim.edgeVertBonds  = new EdgeVertBond[ sim.nEdgeVertBonds ];
    for( int i=0; i<craft.sliders.size(); i++ ){
        Slider* o = craft.sliders[i];
        printf("sliders2edgeverts[%i]: SliderID=%i, ivert=%i, path.n=%i\n", i, o->id, o->ivert, o->path.n);
        if (o->path.n <= 0) { printf("  ERROR: Slider path not initialized! Skipping.\n"); exit(0); }
        if (o->ivert < 0 || o->ivert >= sim.nPoint) { printf("  ERROR: Slider ivert=%i is out of bounds [0..%i]! Skipping.\n", o->ivert, sim.nPoint-1); exit(0); }
        EdgeVertBond& ev = sim.edgeVertBonds[i];
        ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );
        ev.verts.z = o->ivert;
        //ev.K = 1000.0;
        ev.K = 100000.0;
        //ev.K = 0.0;
        o->speed = 1.0;
        //o->speed = 0.1;
        //o->updateEdgeVerts( sim.points );
    }
}

void makeBBoxes( const SpaceCraft& craft, TrussDynamics_d& sim ){
    // // --- Bounding boxes
    int nbuck  = exportBuckets( craft,             0, 16, true );
    //printf( "nbuck %i \n", nbuck );
    sim.recallocBBs( nbuck );
    sim.pointBBs.cleanO2C(0);  // by default we assign all points to cell 0
    int nbuck_ = exportBuckets( craft, &sim.pointBBs, 16, true );
    //sim.pointBBs.printCells();
    //sim.pointBBs.printObjCellMaping(0,100);   // see if obj2cell is incorrect at the beginning
    sim.pointBBs.updateCells();
    //printf("### pointBBs.printCells() \n"); sim.pointBBs.printCells();
    sim.edgesToBBs();
    sim.edgeBBs.updateCells();
    //printf("### edgeBBs.printCells() \n"); sim.edgeBBs.printCells();
    updatePointBBs( sim.pointBBs, sim.BBs, sim.points,            true );  // It crashes here because of the wrong obj2cell mapping
    updateEdgeBBs ( sim.edgeBBs,  sim.BBs, sim.bonds, sim.points, false );
    //updateEdgeBBs ( sim.edgeBBs, sim.BBs, sim.bonds, sim.points, true );
    //sim.printBBs();
}

int makePointCunks( Buckets& ebuck, int2* edges, Buckets& pbuck ){
    //printf( "makePointCunks() ncell %i nobj %i \n", ebuck.ncell, ebuck.nobj );
    std::unordered_set<int> ps; // we use set to remove duplicates
    std::vector<int> p2c;       // we use vector because we do not know the sizes in advance
    pbuck.cellNs = new int[ ebuck.ncell ];
    pbuck.cellI0s= new int[ ebuck.ncell ];
    //DEBUG
    //ebuck.printCells();
    //printf( "makePointCunks() ncell %i nobj %i \n", ebuck.ncell, ebuck.nobj );
    //exit(0);
    // --- select points from edges (without duplicates)
    for(int ib=0; ib<ebuck.ncell; ib++){
        int i0 =      ebuck.cellI0s[ib];
        int i1 = i0 + ebuck.cellNs [ib]; 
        ps.clear();
        //printf( "makePointCunks() ib %i i0 %i i1 %i \n", ib, i0, i1 );
        for(int j=i0; j<i1; j++){
            int ie = ebuck.cell2obj[j];
            int2 e = edges[ie];
            ps.insert( e.x );
            ps.insert( e.y );
        }
        pbuck.cellNs [ib] = ps.size();
        pbuck.cellI0s[ib] = p2c.size();
        for( int p : ps ){ p2c.push_back(p); }  // copy set to vector
    }
    // --- copy point from vector to pbuck
    pbuck.nobj     = p2c.size();
    pbuck.cell2obj = new int[ pbuck.nobj ];
    //pbuck.obj2cel  = new int[ pbuck.nobj ];
    for(int i=0; i<pbuck.nobj; i++){ 
        int j = p2c[i];
        pbuck.cell2obj[i]=j;
        //pbuck.obj2cel [j]=i; 
    }
    return pbuck.nobj;
}

int checkSpaceCraftMesh(const Builder2& mesh, const SpaceCraft& craft, bool bPrint=true, bool bExit=true ){
    std::vector<bool> connected(mesh.verts.size(), false);
    for (const auto& edge : mesh.edges) {
        if (edge.x < mesh.verts.size()) connected[edge.x] = true;
        if (edge.y < mesh.verts.size()) connected[edge.y] = true;
    }
    int unconnected_count = 0;
    for (size_t i = 0; i < mesh.verts.size(); ++i) {
        if (!connected[i]) {
            if (bPrint){ 
                printf("WARNING: checkSpaceCraftMesh() vertex %3i at (%10.3f, %10.3f, %10.3f) is not connected to any edge.\n",  i, mesh.verts[i].pos.x, mesh.verts[i].pos.y, mesh.verts[i].pos.z); 
                printf("        - vertex %i belong to SpaceCraft component: \n", i);
            }
            int j = craft.find_mesh_element(i, 'v', bPrint );
            if(j>=0){ unconnected_count++; }
        }
    }
    // if (unconnected_count > 0) {  
    //     if(bPrint) printf("ERROR: Mesh::Builder2::checkAllPointsConnected() found %i unconnected vertices.\n", unconnected_count);
    //     if(bExit)  exit(0);
    // }
    return unconnected_count;    
}


// Helper function to process Node objects
void node_to_mesh(Node* o, Builder2& mesh) {
    // The bound nodes should not be inserted as new verts, but rather as a reference to existing verts
    if (o->boundTo == 0) {
        o->ivert = mesh.verts.size();
        mesh.vert(o->pos);
    }
}

void nodeBlock_to_mesh(Node* o, Builder2& mesh, SpaceCraft* craft) {
    mesh.block();
    int pre_verts_size = mesh.verts.size();
    int pre_edges_size = mesh.edges.size();
    int pre_chunks_size = mesh.chunks.size();
    if(o->boundTo == 0){
        //Mat3d rot = Mat3dIdentity;
        //printf("nodeBlock_to_mesh() o->boundTo=%i craft=%p\n", o->boundTo, craft );
        printf("nodeBlock_to_mesh() craft->nodeMeshes.size()=%i craft->defaultNodeMesh=%i\n", craft->nodeMeshes.size(), craft->defaultNodeMesh );
        CMesh& nodeMesh = craft->nodeMeshes[ craft->defaultNodeMesh ];
        Vec3d size = Vec3d{1,1,1}*o->size;
        //printf("nodeBlock_to_mesh() nodeMesh.nvert=%i nodeMesh.nedge=%i nodeMesh.nfaces=%i\n", nodeMesh.nvert, nodeMesh.nedge, nodeMesh.nfaces);
        Quat4i i0s = mesh.addCMesh(nodeMesh, true, o->pos, size, 0, o->edge_type );
        //mesh.printChunkRange( i0s.w, mesh.chunks.size() );
    }else{ // attach node using Mesh::Builder2::make_anchor
        mesh.make_anchor_point(o->pos, g_AnchorStickType, g_Rcollapse, g_Ranchor );
    }
    o->pointRange = {pre_verts_size, (int)mesh.verts.size()};
    o->stickRange = {pre_edges_size, (int)mesh.edges.size()};
    o->chunkRange = {pre_chunks_size, (int)mesh.chunks.size()};
}

// Helper function to process Girder objects
void girder_to_mesh(Girder* o, Builder2& mesh) {
    o->update_nodes(); // If girder bound to BoundNode, the vert indexes may need to be updated
    mesh.block();
    mesh.girder1(o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a, o->st);
    mesh.girder1_caps(o->nodes.x->ivert, o->nodes.y->ivert, o->st.x);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

void girderBlocks_to_mesh(Girder* o, Builder2& mesh, SpaceCraft* craft) {
    o->update_nodes(); // If girder bound to BoundNode, the vert indexes may need to be updated
    mesh.block();
    // get chunk ranges
    Node& n1 = *o->nodes.x;
    Node& n2 = *o->nodes.y;
    mesh.bridgeFacingPolygons( n1.pos, n2.pos, n1.chunkRange, n2.chunkRange, o->nseg, o->st );
    //mesh.(o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a, o->st);
    //mesh.girder1_caps(o->nodes.x->ivert, o->nodes.y->ivert, o->st.x);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

// Helper function to process Rope objects
void rope_to_mesh(Rope* o, Builder2& mesh) {
    o->update_nodes();
    mesh.block();
    Quat4i& b = mesh.blocks.back();
    mesh.rope(o->nodes.x->ivert, o->nodes.y->ivert, o->face_mat, o->nseg);
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

// Helper function to process Ring objects
void ring_to_mesh(Ring* o, Builder2& mesh) {
    o->update_nodes();
    mesh.block();
    mesh.wheel(o->pose.pos, o->pose.pos + o->pose.rot.b * o->R, o->pose.rot.c, o->nseg, o->wh, o->st);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

// Helper function to process Slider objects
void slider_to_mesh(Slider* o, Builder2& mesh) {
    // The slider's position 'o->pos' should have been pre-calculated
    // to be the intersection point on its target girder or other component.
    // Record state before creating the anchor point geometry
    int pre_verts_size = mesh.verts.size();
    int pre_edges_size = mesh.edges.size();
    Vec2i host_point_range = o->boundTo->pointRange;
    // Generate the anchor point mesh. This adds a new vertex and connects it to the nearby structure.
    o->ivert = mesh.make_anchor_point(o->pos, g_AnchorStickType, g_Rcollapse, g_Ranchor, 0, 0, host_point_range.x, host_point_range.y);
    // The slider now "owns" the geometry of its anchor point.
    o->pointRange = {pre_verts_size, (int)mesh.verts.size()};
    o->stickRange = {pre_edges_size, (int)mesh.edges.size()};
    printf("slider_to_mesh() o->id %i o->ivert=%i o->pointRange(%i,%i) o->stickRange(%i,%i)\n", o->id, o->ivert, o->pointRange.x, o->pointRange.y, o->stickRange.x, o->stickRange.y);
}

// Helper function to process Radiator objects
void radiator_to_mesh(Radiator* o, Builder2& mesh, SpaceCraft* craft) {
    mesh.block();
    o->print(); // Assuming print() is for debugging and can stay
    const Girder& g1 = *craft->girders[o->g1];
    const Girder& g2 = *craft->girders[o->g2];
    mesh.plateOnGriders({10, 1}, g1.pointRange, g2.pointRange, {4, 4}, {-1, -1}, o->g1span, o->g2span, {0, 1, 2, 3});
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

// Helper function to process Weld objects
void weld_to_mesh(Weld* o, Builder2& mesh) {
    mesh.block();
    mesh.bondsBetweenVertRanges(o->comps.x->pointRange, o->comps.y->pointRange, o->Rmax, o->face_mat);
    Quat4i& b = mesh.blocks.back();
    o->pointRange = {b.x, (int)mesh.verts.size()};
    o->stickRange = {b.y, (int)mesh.edges.size()};
}

void BuildCraft_truss( Builder2& mesh, SpaceCraft& craft, double max_size=-1 ){
    //printf( "BuildCraft_truss() \n" );
    if(max_size>0){ mesh.max_size=max_size; };
    int i=0;
    mesh.block();
    int ip0 = mesh.verts.size();
    //printf("BuildCraft_truss().nodes n=%i\n", craft.nodes.size() );
    for(Node* o: craft.nodes){ node_to_mesh(o, mesh);  }
    //printf("BuildCraft_truss().girders n=%i\n", craft.girders.size() );
    for(Girder* o: craft.girders){ girder_to_mesh(o, mesh); i++; }
    //printf("BuildCraft_truss().ropes n=%i\n", craft.ropes.size() );
    for(Rope* o: craft.ropes){ rope_to_mesh(o, mesh);}
    // --- Rings
    //printf("BuildCraft_truss().rings n=%i\n", craft.rings.size() );
    for(Ring* o: craft.rings){ ring_to_mesh(o, mesh);}
    // --- Radiators
    //printf("BuildCraft_truss().radiators n=%i\n", craft.radiators.size() );
    for(Radiator* o : craft.radiators ){ radiator_to_mesh(o, mesh, &craft); };
    //printf("BuildCraft_truss().Welds n=%i\n", craft.welds.size() );
    for( Weld* o : craft.welds ){ weld_to_mesh(o, mesh);};
    // --- Shields
    //for( const Shield* o : craft.shields ){
    //    //mesh.block();
    //    //drawPlate_mesh(mesh, o, nodes.data(), girders.data() );
    //};
    //int plate_quad( int ip00, int ip01, int ip10, int ip11, Quat4i typs={-1,-1,-1,-1}, Vec2i n={1,1}, int fillType=1 );
    // // ToDo: Shields
    // // ToDo: Radiators
    // // ToDo: Thrusters
    // // ToDo: Tanks
    // // ToDo: maybe Tanks would be better simulated as rigid-body objects, but than we need to couple the Truss (SoftBody) and RigidBody ?
    // if( bTanks ){
    //     for(Tank o: tanks){
    //         //printf("DEBUG toTruss : tank #%i \n", i);
    //         Vec3d p0 = o.pose.pos + o.pose.rot.a*o.span.a*0.5;
    //         Vec3d p1 = o.pose.pos - o.pose.rot.a*o.span.a*0.5;
    //         int edgeMat = workshop->panelMaterials.vec[ o.face_mat ]->stickMaterialId;
    //         truss.makeCylinder( p0, p1, o.span.b, o.span.b, -1, -1, 1.0, o.face_mat, o.face_mat );
    //         Vec2i& bak = truss.blocks.back();
    //         o.pointRange  = {bak.x,truss.points.size()};
    //         o.stickRange = {bak.y,truss.edges .size()};
    //         i++;
    //     }
    // }
    printf( "BuildCraft_truss() DONE! : npoint %i nstick %i nblocks %i \n", mesh.verts.size(), mesh.edges.size(), mesh.blocks.size()  );
}

// New block-based construction function
void BuildCraft_blocks( Builder2& mesh, SpaceCraft& craft, double max_size=-1, double node_scale=1.0 ){
    printf( "### BuildCraft_blocks()\n" );
    if(max_size>0){ mesh.max_size=max_size; };
    // ------ PASS 1: Generate Node Blocks (Octahedra) -------
    mesh.block();
    int ip0 = mesh.verts.size();
    printf("BuildCraft_blocks().nodes n=%i\n", craft.nodes.size());
    for(Node* o: craft.nodes){ 
        if (o->component_kind() == (int)ComponetKind::Slider){ 
            //printf("BuildCraft_blocks().nodes: Slider processed in BuildCraft_blocks().nodes loop \n");
            //exit(0);
            continue;
        }
        nodeBlock_to_mesh(o, mesh, &craft); 
    }
    printf("BuildCraft_blocks().girders n=%i\n", craft.girders.size() );   
    for(Girder* o: craft.girders){ girderBlocks_to_mesh(o, mesh, &craft); }
    // ------ PASS 3: Generate Rings -------
    printf("BuildCraft_blocks().rings n=%i\n", craft.rings.size() );
    for(Ring* o: craft.rings){ ring_to_mesh(o, mesh); }
    // ------ PASS 4: Generate Slider Anchors -------
    printf("BuildCraft_blocks().sliders n=%i\n", craft.sliders.size() );
    for(Slider* o: craft.sliders){ slider_to_mesh(o, mesh); }
    // ------ PASS 5: Generate Ropes -------
    printf("BuildCraft_blocks().ropes n=%i\n", craft.ropes.size() );
    for(Rope* o: craft.ropes){ rope_to_mesh(o, mesh);}

    // --- ToDo Future ---
    //printf("BuildCraft_blocks().radiators n=%i\n", craft.radiators.size() );
    
    //for(Radiator* o : craft.radiators ){ radiator_to_mesh(o, mesh, &craft); };
    //for( Weld* o : craft.welds ){ weld_to_mesh(o, mesh);};
    
    // --- Final integrity check
    //mesh.checkAllPointsConnected(craft, true, true);
    checkSpaceCraftMesh(mesh, craft);

    printf("BuildCraft_blocks() DONE! : npoint %i nstick %i nblocks %i\n", mesh.verts.size(), mesh.edges.size(), mesh.blocks.size());
}

void makeTrussShape( Mesh::Builder2& mesh, int ishape, int nseg, double R, double r, Mat3d rot=Mat3dIdentity, Vec3d p0=Vec3dZero, double Rcolapse=0.1, double Ranchor=10.0 ){
    printf("makeTrussShape(%i,%i,%f,%f)\n",ishape,nseg,R,r);
    //StickMaterial *o = new StickMaterial();
    //shop.add_Material     ( "Steel", 7.89e+3, 1.2e+9, 1.2e+9, 200.0e+9, 200.0e+9, 0.85, 800 );
    //shop.add_StickMaterial( "GS1_long", "Steel", 0.1, 0.005, 0.0 );
    Quat4i stickTypes{0,0,0,0};
    Vec3d p1 = rot.a * R;
    Vec3d up = rot.b;
    Vec3d ax = rot.c;
    mesh.clear();
    switch(ishape){
        case 0: mesh.ngon          ( p0, p1                            , ax, nseg,             stickTypes.x       ); break;
        case 1: mesh.wheel         ( p0, p1                            , ax, nseg, Vec2d{r,r}, stickTypes         ); break;
        case 2: mesh.girder1       ( p0, (p1-p0).normalized()*r*4*nseg , ax, nseg, r,          stickTypes,   true ); break;
        case 3: mesh.triangle_strip( p0, (p1-p0).normalized()*r*nseg   , up, nseg, r,          stickTypes.x, true ); break;
        case 4: mesh.rope          ( p0, (p1-p0).normalized()*r*nseg   ,     nseg,             stickTypes.x, stickTypes.y, Rcolapse, Ranchor ); break;
    }
    mesh.printSizes();
    //if(bDouble){ to_TrussDynamics ( sim2,   mesh,        p0,         ax      ); }
    //if(bSingle){ to_OCL_Orb( sim_cl, mesh, (Vec3f)p0, (Vec3f)(ax*5.0) ); }
    //distort_points( sim_cl.nPoint, sim_cl.points, sim_cl.points, 1.0, Vec3fOne, 15454 ); 
};

} // namespace SpaceCrafting

#endif
