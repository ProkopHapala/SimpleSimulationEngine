
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

using namespace SpaceCrafting;

namespace Mesh{

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
        if(e.w>=shop.stickMaterials.vec.size()){ printf( "ERROR in exportSim() mesh.edges[%i].type=%i > stickMaterials.size()\n", i, e.w, e.w>=shop.stickMaterials.vec.size() ); exit(0); }
        const StickMaterial& mat = *shop.stickMaterials.vec[e.w];
        stickTypes.insert(e.w);
        // l0, kPress, kPull, damping
        //double l0 = (mesh.verts[e.y].pos - mesh.verts[e.x].pos ).norm();
        double l0 = (sim.points[e.y].f - sim.points[e.x].f ).norm() * ( 1.f - mat.preStrain );
        double mass = l0*mat.linearDensity;
        sim.points[e.x].w += mass*0.5;
        sim.points[e.y].w += mass*0.5;
        Quat4d param = (Quat4d){ l0, mat.Kpush/l0, mat.Kpull/l0, mat.damping }; 


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
    //printf( "reloadShip().updateSlidersPaths \n" );
    // update ring slider paths
    for( Ring* o : craft.rings ){
        o->updateSlidersPaths( true, true, sim.points );
    }
    //printf( "reloadShip().SlidersToEdgeVerts \n" );
    // ----- Sliders to EdgeVerts
    sim.nEdgeVertBonds = craft.sliders.size();
    sim.edgeVertBonds  = new EdgeVertBond[ sim.nEdgeVertBonds ];
    for( int i=0; i<craft.sliders.size(); i++ ){
        Slider* o = craft.sliders[i];
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


void BuildCraft_truss( Builder2& mesh, SpaceCraft& craft, double max_size=-1 ){
    //printf( "BuildCraft_truss() \n" );
    if(max_size>0){ mesh.max_size=max_size; };
    int i=0;
    mesh.block();
    int ip0 = mesh.verts.size();
    for(Node* o: craft.nodes){
        // the bound nodes should not be inserted as a new verts, but rather as a reference to existing verts
        if( o->boundTo == 0 ){ 
            o->ivert = mesh.verts.size();
            mesh.vert( o->pos );
        }
    }
    //printf("BuildCraft_truss().girders n=%i\n", craft.girders.size() );
    for(Girder* o: craft.girders){
        //printf("DEBUG toTruss : girder #%i \n", i);
        o->update_nodes(); // if girder bound to BoundNode, the vert indexes may need to be updated
        mesh.block();
        mesh.girder1(  o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a, o->st );
        mesh.girder1_caps( o->nodes.x->ivert, o->nodes.y->ivert, o->st.x );
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //o.print();
        //printf( "BuildCraft_truss() girder.pointRange(%i,%i)\n", o.pointRange.x, o.pointRange.y );
        i++;
    }
    //printf("BuildCraft_truss().ropes n=%i\n", craft.ropes.size() );
    for(Rope* o: craft.ropes){
        o->update_nodes();
        // Node* nd1 = o->nodes.x;
        // Node* nd2 = o->nodes.y;
        // if( (nd1->boundTo!=0)||(nd2->boundTo!=0) ){ 
        //     printf("Bound_Rope[%i] to mesh\n", o->id );
        //     printf("  node.x ivert=%i pos(%g,%g,%g) boundTo: \n", nd1->ivert, nd1->pos.x,nd1->pos.y,nd1->pos.z ); if(nd1->boundTo)nd1->boundTo->print();
        //     printf("  node.y ivert=%i pos(%g,%g,%g) boundTo: \n", nd2->ivert, nd2->pos.x,nd2->pos.y,nd2->pos.z ); if(nd2->boundTo)nd2->boundTo->print();
        // }
        mesh.block(); Quat4i& b = mesh.blocks.back();
        //truss.edges.push_back( (TrussEdge){o.p0,o.p1,0} );
        mesh.rope( o->nodes.x->ivert,o->nodes.y->ivert, o->face_mat, o->nseg );
        
        // // Rope Dampers ( to dampe perpendicular oscillations )
        // double damper_length = 10.0;
        // Vec3d d = o->nodes.y->pos - o->nodes.x->pos; d.normalize();
        // Vec3d up,lf; d.getSomeOrtho( up, lf );
        // up.mul( damper_length );
        // for(int i=0; i<o->nseg-1; i++){
        //     int ov = b.x+i;
        //     int v  = mesh.vert( mesh.verts[ ov ].pos + up ); // duplicate the vert
        //     mesh.edge( ov,v, o->face_mat );
        // }
        
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
    }
    // --- Rings
    //printf("BuildCraft_truss().rings n=%i\n", craft.rings.size() );
    i=0;
    for(Ring* o: craft.rings){
        o->update_nodes();

        // Node** nd = (Node**)&o->nodes;
        // printf("Ring[%i] to mesh\n", o->id );
        // for(int i=0; i<4; i++){
        //     Node* n = nd[i]; printf("  node[%i] ivert=%i pos(%g,%g,%g) boundTo: \n", i, n->ivert, n->pos.x,n->pos.y,n->pos.z ); if(n->boundTo)n->boundTo->print();
        // }

        mesh.block();
        //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
        mesh.wheel( o->pose.pos, o->pose.pos+o->pose.rot.b*o->R, o->pose.rot.c, o->nseg, o->wh, o->st );
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //printf( "BuildCraft_truss() ring.pointRange(%i,%i)\n", o.pointRange.x, o.pointRange.y );
        i++;
    }
    // --- Radiators
    //printf("BuildCraft_truss().radiators n=%i\n", craft.radiators.size() );
    for(Radiator* o : craft.radiators ){
        mesh.block();
        o->print();
        const Girder& g1 = *craft.girders[o->g1];
        const Girder& g2 = *craft.girders[o->g2];
        mesh.plateOnGriders( {10,1}, g1.pointRange, g2.pointRange, {4,4}, {-1,-1}, o->g1span, o->g2span, {0,1,2,3} );
        //break;
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //break;
    };
    //printf("BuildCraft_truss().Welds n=%i\n", craft.welds.size() );
    for( Weld* o : craft.welds ){
        //printf("Welds[]\n");
        //o->print();
        mesh.block();
        mesh.bondsBetweenVertRanges( o->comps.x->pointRange, o->comps.y->pointRange, o->Rmax, o->face_mat );
        //break;
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //break;
    };
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

void makeTrussShape( Mesh::Builder2& mesh, int ishape, int nseg, double R, double r, Mat3d rot=Mat3dIdentity, Vec3d p0=Vec3dZero ){
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
        case 4: mesh.rope          ( p0, (p1-p0).normalized()*r*nseg   ,     nseg,             stickTypes.x       ); break;
    }
    mesh.printSizes();
    //if(bDouble){ to_TrussDynamics ( sim2,   mesh,        p0,         ax      ); }
    //if(bSingle){ to_OCL_Orb( sim_cl, mesh, (Vec3f)p0, (Vec3f)(ax*5.0) ); }
    //distort_points( sim_cl.nPoint, sim_cl.points, sim_cl.points, 1.0, Vec3fOne, 15454 ); 
};

} // namespace SpaceCrafting

#endif
