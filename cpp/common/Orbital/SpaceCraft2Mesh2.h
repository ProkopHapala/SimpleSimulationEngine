
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

using namespace SpaceCrafting;

namespace Mesh{

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
    printf( "Mesh::girder1() n=%i p0(%g,%g,%g) p1(%g,%g,%g) up(%g,%g,%g) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, up.x,up.y,up.z );
    Vec3d dir = p1-p0;
    double length = dir.normalize();
    up.makeOrthoU(dir);
    Vec3d side; side.set_cross(dir,up);
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
int wheel( Builder2& mesh, Vec3d p0, Vec3d p1, Vec3d ax, int n, double width, Quat4i stickTypes ){
    printf( "Truss::wheel() n=%i p0(%g,%g,%g) p1(%g,%g,%g) ax(%g,%g,%g) \n", n, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, ax.x,ax.y,ax.z );
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
        mesh.vert( p0 + R*r );
        //mesh.vert( p0 + R*(r+width) );
        mesh.vert( p0 + R*(r-width) );
        // points.push_back( p0 +  R*(r-width) );
        // points.push_back( p0 +  R*(r+width) );
        rot.mul_cmplx(drot);
        R       = dir*rot.a + side*rot.b;
        mesh.vert( p0 + ax*-width + R*r );
        mesh.vert( p0 + ax*-width + R*r );
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
    printf( "Mesh::panel() n(%i,%i) w=%g p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g) \n", n.x,n.y, p00.x,p00.y,p00.z, p01.x,p01.y,p01.z, p10.x,p10.y,p10.z, p11.x,p11.y,p11.z );
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

void BuildCraft_truss( Builder2& mesh, const SpaceCraft& craft ){
    printf( "BuildCraft_truss() \n" );
    int i=0;
    mesh.block();
    int ip0 = mesh.verts.size();
    for(Node o: craft.nodes){
        mesh.vert( o.pos );
    }
    for(Rope o: craft.ropes){
        //truss.edges.push_back( (TrussEdge){o.p0,o.p1,0} );
        mesh.rope( o.p0,o.p1, o.face_mat );
    }
    for(Girder o: craft.girders){
        //printf("DEBUG toTruss : girder #%i \n", i);
        girder1( mesh, craft.nodes[o.p0].pos, craft.nodes[o.p1].pos, o.up, o.nseg, o.wh.a, o.st );
        girder1_caps( mesh, o.p0, o.p1, o.st.x );
        Quat4i& b = mesh.blocks.back();
        o.poitRange  = {b.x,(int)mesh.verts.size()};
        o.stickRange = {b.y,(int)mesh.edges.size()};
        i++;
    }
    i=0;
    for(Ring o: craft.rings){
        //printf("DEBUG toTruss : ring #%i  %f   %f \n", i, o.nseg, o.wh.a );
        wheel( mesh, o.pose.pos, o.pose.pos+o.pose.rot.b*o.R, o.pose.rot.c, o.nseg, o.wh.a, o.st );
        Quat4i& b = mesh.blocks.back();
        o.poitRange  = {b.x,(int)mesh.verts.size()};
        o.stickRange = {b.y,(int)mesh.edges.size()};
        i++;
    }
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
            o.poitRange  = {bak.x,truss.points.size()};
            o.stickRange = {bak.y,truss.edges .size()};
            i++;
        }
    }
    */
    printf( "BuildCraft_truss() DONE! : npoint %i nstick %i nblocks %i \n", mesh.verts.size(), mesh.edges.size(), mesh.blocks.size()  );
}


} // namespace SpaceCrafting

#endif
