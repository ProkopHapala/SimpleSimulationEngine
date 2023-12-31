
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
#include "OrbSim.h"

using namespace SpaceCrafting;

namespace Mesh{

void exportSim( OrbSim_f& sim, const Builder2& mesh, const SpaceCraftWorkshop& shop ){
    {
        printf("#====== Materials \n");      for(const Material*      o: shop.materials.vec     ){  o->print();    }
        //printf("#====== PanelMaterials \n"); for(const PanelMaterial* o: shop.panelMaterials.vec){  o->print();    }
        printf("#====== StickMaterials \n"); for(const StickMaterial* o: shop.stickMaterials.vec){  o->print();  }
        //exit(0);
    }
    printf( "exportSim() START \n" );
    int np = mesh.verts.size();
    int nb = mesh.edges.size();
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
    for(int i=0; i<=nneighmax; i++){ printf( "exportSim() nneigh_hist[%i] %i \n", i, nneigh_hist[i] ); }
    //exit(0);

    //printf( "exportSim() nneighmax %i \n", nneighmax );
    sim.recalloc(np, nneighmax, nb );
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
        double l0 = (sim.points[e.y].f - sim.points[e.x].f ).norm();
        double mass = l0*mat.linearDensity;
        sim.points[e.x].w += mass*0.5;
        sim.points[e.y].w += mass*0.5;
        Quat4f param = (Quat4f){ l0, mat.Kpush/l0, mat.Kpull/l0, mat.damping }; 
        {
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
    printf( "exportSim() DONE! \n" );
    delete [] nneighs;
    //exit(0);
}

/**
 * Creates a girder in the truss structure.
 * 
 * @param p0 starting point
 * @param p1 ending point
 * @param up up-vector 
 * @param n  number of segments along the girder 
 * @param width The width of the girder.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int girder1( Builder2& mesh, Vec3d p0, Vec3d p1, Vec3d up, int n, double width, Quat4i stickTypes ){
    // ToDo: ad p0,p1 etc. - maybe we should rather specify indexes of existing verts rather than positions of new verts ?
    //                       No. this is done in girder1_caps
    //printf( "Mesh::girder1() n=%i p0(%g,%g,%g) p1(%g,%g,%g) up(%g,%g,%g) stickTypes(%i,%i,%i,%i)\n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, up.x,up.y,up.z,   stickTypes.x,stickTypes.y,stickTypes.z,stickTypes.w );
    Vec3d dir = p1-p0;
    double length = dir.normalize();
    up.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,up); side.normalize();
    //print(dir); print(up); print(side); printf("dir up side \n");
    double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    int dnp = 4;
    int i00   = mesh.verts.size();
    //int ibloc = mesh.block();   // this is better to call manually from outside
    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;
        mesh.vert( p0 + side*-width + dir*(dl*(1+2*i  )) );
        mesh.vert( p0 + side*+width + dir*(dl*(1+2*i  )) );
        mesh.vert( p0 + up  *-width + dir*(dl*(1+2*i+1)) );
        mesh.vert( p0 + up  *+width + dir*(dl*(1+2*i+1)) );
        mesh.edge( i00,i01,stickTypes.y );
        mesh.edge( i10,i11,stickTypes.y );
        mesh.edge( i00,i10,stickTypes.z );
        mesh.edge( i00,i11,stickTypes.z );
        mesh.edge( i01,i10,stickTypes.z );
        mesh.edge( i01,i11,stickTypes.z );
        if( i<(n-1) ){
            mesh.edge( i10,i00+dnp,stickTypes.w );
            mesh.edge( i10,i01+dnp,stickTypes.w );
            mesh.edge( i11,i00+dnp,stickTypes.w );
            mesh.edge( i11,i01+dnp,stickTypes.w );
            mesh.edge( i00,i00+dnp,stickTypes.x );
            mesh.edge( i01,i01+dnp,stickTypes.x );
            mesh.edge( i10,i10+dnp,stickTypes.x );
            mesh.edge( i11,i11+dnp,stickTypes.x );
        }
        i00+=dnp;
    }
    //return ibloc;
    return i00;
}



/* @brief : plot planar fill between two strips of vertexes (e.g. between two girders or ropes)
*
*/
// int plateOnGriders( Builder2& mesh, Vec2i ns, Vec2i iblocks, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes ){
int plateOnGriders( Builder2& mesh, Vec2i ns, Vec2i prange1, Vec2i prange2, Vec2i byN, Vec2i offs, Vec2d span1, Vec2d span2, Quat4i stickTypes ){
    //printf( "plateOnGriders() ns(%i,%i) iblocks(%i,%i) byN(%i,%i) offs(%i,%i) span1(%6.3f,%6.3f) span2(%6.3f,%6.3f) \n", ns.x,ns.y, iblocks.x,iblocks.y, byN.x,byN.y,  offs.x,offs.y, span1.x,span1.y,  span2.x,span2.y );
    int n1   = (prange1.y - prange1.x)/byN.x;
    int n2   = (prange2.y - prange2.x)/byN.y;
    double step = 1.0/(ns.x-1);
    if(offs.x<0){
        int i0 = prange1.x + byN.x*(n1/2);
        offs.x=mesh.findClosestVert( mesh.verts[ prange2.x + byN.y*(n2/2) ].pos,i0,byN.x);
    }
    if(offs.y<0){
        int i0 = prange2.x + byN.y*(n2/2);
        offs.y=mesh.findClosestVert( mesh.verts[ prange1.x + byN.x*(n1/2) ].pos,i0,byN.y);
    }
    for(int i=0;i<ns.x; i++){
        double c = i*step;
        int i1 = prange1.x + offs.x + byN.x*(int)(  ( span1.x*(1-c)+span1.y*c)* n1 + 0.5 );
        int i2 = prange2.x + offs.y + byN.y*(int)(  ( span2.x*(1-c)+span2.y*c)* n2 + 0.5 );
        //printf( "plateOnGriders()[%i] (%i/%i) (%i/%i)\n", i, i1,n1, i2,n2 );
        mesh.edge( i1,i2,stickTypes.x );
    }
    //return ibloc;
    return 0;
}




int girder1_caps( Builder2& mesh, int ip0, int ip1, int kind ){
    //printf( "Truss::girder1_caps() ip0=%i ip1=%i ipbeg=%i ipend=%i \n", ip0, ip1, ipbeg, ipend );
    // it is expected that we call this immediately after girder1 - so we know the indexes of the first and last point of the girder block
    int ipbeg = mesh.blocks.back().x;  // index of the first point of the girder block
    int ipend = mesh.verts.size()-4;   // index of the last point  of the girder block
    int ie0 = mesh.edges.size();
    mesh.edge( ip0,ipbeg+0,kind );
    mesh.edge( ip0,ipbeg+1,kind );
    mesh.edge( ip0,ipbeg+2,kind );
    mesh.edge( ip0,ipbeg+3,kind );
    mesh.edge( ip1,ipend+0,kind );
    mesh.edge( ip1,ipend+1,kind );
    mesh.edge( ip1,ipend+2,kind );
    mesh.edge( ip1,ipend+3,kind );
    return ie0;
}

int girder1( Builder2& mesh, int ip0, int ip1, Vec3d up, int n, double width, Quat4i stickTypes ){
    int ib = girder1( mesh, mesh.verts[ip0].pos, mesh.verts[ip1].pos, up, n, width, stickTypes );
    girder1_caps( mesh, ip0, ip1, stickTypes.x );
    return ib;
}

/**
 * Creates a wheel-shaped truss structure.
 *
 * @param p0 starting point ( center of the wheel )
 * @param p1 ending point   ( on the perimeter of the wheel )
 * @param ax axis vector
 * @param n The number of segments along the perimeter of the wheel.
 * @param width The width of the wheel rim.
 * @return The index of the created block containing the starting indexes of the points and edges
 */
int wheel( Builder2& mesh, Vec3d p0, Vec3d p1, Vec3d ax, int n, Vec2d wh, Quat4i stickTypes ){
    //printf( "Truss::wheel() n=%i p0(%g,%g,%g) p1(%g,%g,%g) ax(%g,%g,%g) stickTypes(%i,%i,%i,%i) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, ax.x,ax.y,ax.z,   stickTypes.x,stickTypes.y,stickTypes.z,stickTypes.w );
    //int kind_long   = 0;
    //int kind_perp   = 1;
    //int kind_zigIn  = 2;
    //int kind_zigOut = 3;
    Vec3d dir = p1-p0;
    double r = dir.normalize();
    ax.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,ax);
    //double dl = length/(2*n + 1);
    //int dnb = 2+4+4+4;
    //print(dir); print(ax); print(side); printf("dir up side \n");
    int dnp = 4;
    int i00 = mesh.verts.size();
    int i000 = i00;

    Vec2d  rot = {1.0,0.0};
    Vec2d drot; drot.fromAngle( M_PI/n );
    //int ibloc = mesh.blocks.size();
    //mesh.blocks.push_back( {i00,edges.size()} );
    for (int i=0; i<n; i++){
        int i01=i00+1; int i10=i00+2; int i11=i00+3;

        Vec3d R = dir*rot.a + side*rot.b;
        //mesh.vert( p0 + R*r );
        mesh.vert( p0 + R*(r+wh.x) );
        mesh.vert( p0 + R*(r-wh.x) );
        // points.push_back( p0 +  R*(r-width) );
        // points.push_back( p0 +  R*(r+width) );
        rot.mul_cmplx(drot);
        R       = dir*rot.a + side*rot.b;
        mesh.vert( p0 + ax*+wh.y + R*r );
        mesh.vert( p0 + ax*-wh.y + R*r );
        // points.push_back( p0 + ax*-width + R*r );
        // points.push_back( p0 + ax*+width + R*r );
        rot.mul_cmplx(drot);

        mesh.edge( i00,i01,stickTypes.y );
        mesh.edge( i10,i11,stickTypes.y );
        mesh.edge( i00,i10,stickTypes.z );
        mesh.edge( i00,i11,stickTypes.z );
        mesh.edge( i01,i10,stickTypes.z );
        mesh.edge( i01,i11,stickTypes.z );
        // edges .push_back( (TrussEdge){i00,i01,kind_perp}  );
        // edges .push_back( (TrussEdge){i10,i11,kind_perp}  );
        // edges .push_back( (TrussEdge){i00,i10,kind_zigIn} );
        // edges .push_back( (TrussEdge){i00,i11,kind_zigIn} );
        // edges .push_back( (TrussEdge){i01,i10,kind_zigIn} );
        // edges .push_back( (TrussEdge){i01,i11,kind_zigIn} );
        if( i<(n-1) ){
            mesh.edge( i10,i00+dnp,stickTypes.w );
            mesh.edge( i10,i01+dnp,stickTypes.w );
            mesh.edge( i11,i00+dnp,stickTypes.w );
            mesh.edge( i11,i01+dnp,stickTypes.w );
            mesh.edge( i00,i00+dnp,stickTypes.x );
            mesh.edge( i01,i01+dnp,stickTypes.x );
            mesh.edge( i10,i10+dnp,stickTypes.x );
            mesh.edge( i11,i11+dnp,stickTypes.x );
            // edges.push_back( (TrussEdge){i10,i00+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i10,i01+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i00+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i01+dnp,kind_zigOut} );
            // edges.push_back( (TrussEdge){i00,i00+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i01,i01+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i10,i10+dnp,kind_long} );
            // edges.push_back( (TrussEdge){i11,i11+dnp,kind_long} );
        }else{
            mesh.edge( i10,i000+0,stickTypes.w );
            mesh.edge( i10,i000+1,stickTypes.w );
            mesh.edge( i11,i000+0,stickTypes.w );
            mesh.edge( i11,i000+1,stickTypes.w );
            mesh.edge( i00,i000+0,stickTypes.x );
            mesh.edge( i01,i000+1,stickTypes.x );
            mesh.edge( i10,i000+2,stickTypes.x );
            mesh.edge( i11,i000+3,stickTypes.x );
            // edges.push_back( (TrussEdge){i10,i000+0,kind_zigOut} );
            // edges.push_back( (TrussEdge){i10,i000+1,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i000+0,kind_zigOut} );
            // edges.push_back( (TrussEdge){i11,i000+1,kind_zigOut} );
            // edges.push_back( (TrussEdge){i00,i000+0,kind_long} );
            // edges.push_back( (TrussEdge){i01,i000+1,kind_long} );
            // edges.push_back( (TrussEdge){i10,i000+2,kind_long} );
            // edges.push_back( (TrussEdge){i11,i000+3,kind_long} );
        }
        i00+=dnp;
    }
    return i000;
}

/**
 * Creates a panel of truss elements between four corner points. Adds the points and edges to the Truss object.
 * // TODO: make also triangular panel
 * 
 * @param p00 The first corner point.
 * @param p01 The second corner point.
 * @param p10 The third corner point.
 * @param p11 The fourth corner point.
 * @param n The number of subdivisions along each side of the panel.
 * @param width The width of the truss elements.
 */
int panel( Builder2& mesh, Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width, Quat4i stickTypes ){
    // ToDo: ad p00,p01,p10,p11 etc. - maybe we should rather specify indexes of existing verts rather than positions of new verts ?
    //printf( "Mesh::panel() n(%i,%i) w=%g p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g) \n", n.x,n.y, p00.x,p00.y,p00.z, p01.x,p01.y,p01.z, p10.x,p10.y,p10.z, p11.x,p11.y,p11.z );
    //int kind_long   = 0;
    //int kind_perp   = 1;
    //int kind_zigIn  = 2;
    //int kind_zigOut = 3;
    Vec2d step = {1.0/n.a,1.0/n.b};
    int di = 2*n.a-1;
    //int ibloc = mesh.block();   // this is better to call manually from outside
    int i0  = mesh.verts.size();
    //int i00 = mesh.verts.size();
    for (int ib=0; ib<n.b; ib++){
        double db,mb;
        db = ib*step.b;     mb=1-db;
        Vec3d p0  = p00*mb + p10*db;
        Vec3d p1  = p01*mb + p11*db;
        db += 0.5*step.b; mb=1-db;
        Vec3d p0_ = p00*mb + p10*db;
        Vec3d p1_ = p01*mb + p11*db;
        for (int ia=0; ia<n.a; ia++){
            double da,ma;
            da = ia*step.a; ma = 1-da;
            Vec3d p   = p0 *ma + p1 *da;
            //points.push_back( p              );
            mesh.vert( p );
            int bi = i0+di; if( ib==n.b-2 )bi-=ia;
            int dia = 2;    if( ib==n.b-1 )dia=1;
            if (ia<(n.a-1)) mesh.edge( i0,i0+dia,stickTypes.y );
            if (ib<(n.b-1)) mesh.edge( i0,bi    ,stickTypes.y );
            if( (ia<(n.a-1))&&(ib<(n.b-1)) ){ // diagonal
                Vec3d p_  = p0_*ma + p1_*da;
                da += 0.5*step.a; ma=1-da;
                Vec3d p__ = p0_*ma + p1_*da;
                Vec3d up; up.set_cross( p_-p, p__-p ); up.normalize();
                //points.push_back( p__ + up*width );
                mesh.vert( p__ + up*width );
                if( ia<(n.a-2) ) mesh.edge( i0+1,i0+1+dia,stickTypes.z );
                if( ib<(n.b-2) ) mesh.edge( i0+1,bi+1    ,stickTypes.z );
                mesh.edge( i0+1,i0     ,stickTypes.w );
                mesh.edge( i0+1,i0+dia ,stickTypes.w );
                mesh.edge( i0+1,bi     ,stickTypes.w );
                if( ib==n.b-2 )dia=1;
                mesh.edge( i0+1,bi+dia ,stickTypes.w );
                i0++;
            }
            i0++;
        }
    }
    return i0;
}

void BuildCraft_truss( Builder2& mesh, SpaceCraft& craft, double max_size=-1 ){
    printf( "BuildCraft_truss() \n" );
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
    for(Rope* o: craft.ropes){
        o->update_nodes();
        mesh.block();
        //truss.edges.push_back( (TrussEdge){o.p0,o.p1,0} );
        mesh.rope( o->nodes.x->id,o->nodes.y->id, o->face_mat, o->nseg );
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
    }
    for(Girder* o: craft.girders){
        //printf("DEBUG toTruss : girder #%i \n", i);
        o->update_nodes(); // if girder bound to BoundNode, the vert indexes may need to be updated
        mesh.block();
        girder1( mesh, o->nodes.x->pos, o->nodes.y->pos, o->up, o->nseg, o->wh.a, o->st );
        girder1_caps( mesh, o->nodes.x->ivert, o->nodes.y->ivert, o->st.x );
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //o.print();
        //printf( "BuildCraft_truss() girder.pointRange(%i,%i)\n", o.pointRange.x, o.pointRange.y );
        i++;
    }
    // --- Rings
    i=0;
    for(Ring* o: craft.rings){
        o->update_nodes();
        mesh.block();
        //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
        wheel( mesh, o->pose.pos, o->pose.pos+o->pose.rot.b*o->R, o->pose.rot.c, o->nseg, o->wh, o->st );
        Quat4i& b = mesh.blocks.back();
        o->pointRange = {b.x,(int)mesh.verts.size()};
        o->stickRange = {b.y,(int)mesh.edges.size()};
        //printf( "BuildCraft_truss() ring.pointRange(%i,%i)\n", o.pointRange.x, o.pointRange.y );
        i++;
    }
    // --- Radiators
    printf("BuildCraft_truss().radiators\n");
    for(Radiator* o : craft.radiators ){
        mesh.block();
        o->print();
        const Girder& g1 = *craft.girders[o->g1];
        const Girder& g2 = *craft.girders[o->g2];
        plateOnGriders( mesh, {10,1}, g1.pointRange, g2.pointRange, {4,4}, {-1,-1}, o->g1span, o->g2span, {0,1,2,3} );
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
    /*
    // ToDo: Shields
    // ToDo: Radiators
    // ToDo: Thrusters
    // ToDo: Tanks
    // ToDo: maybe Tanks would be better simulated as rigid-body objects, but than we need to couple the Truss (SoftBody) and RigidBody ?
    if( bTanks ){
        for(Tank o: tanks){
            //printf("DEBUG toTruss : tank #%i \n", i);
            Vec3d p0 = o.pose.pos + o.pose.rot.a*o.span.a*0.5;
            Vec3d p1 = o.pose.pos - o.pose.rot.a*o.span.a*0.5;
            int edgeMat = workshop->panelMaterials.vec[ o.face_mat ]->stickMaterialId;
            truss.makeCylinder( p0, p1, o.span.b, o.span.b, -1, -1, 1.0, o.face_mat, o.face_mat );
            Vec2i& bak = truss.blocks.back();
            o.pointRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
    }
    */
    printf( "BuildCraft_truss() DONE! : npoint %i nstick %i nblocks %i \n", mesh.verts.size(), mesh.edges.size(), mesh.blocks.size()  );
}


} // namespace SpaceCrafting

#endif
