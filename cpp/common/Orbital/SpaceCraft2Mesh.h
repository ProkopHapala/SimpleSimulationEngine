
#ifndef  SpaceCraft2Mesh_h
#define  SpaceCraft2Mesh_h

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//#include "Noise.h"
//#include "SphereSampling.h"
//#include "DrawSphereMap.h"

#include "Truss.h"
#include "SpaceCraft.h"
//#include "EditSpaceCraft.h"
#include "MeshBuilder.h"
#include "MeshBuilder2.h"
#include "MeshBuilderDraw.h"

namespace SpaceCrafting{

/**
 * Draws a truss into MeshBuilder object. You can specify the width of the truss tubes. If the width is negative, the truss will be drawn as lines.
 * 
 * @param mesh   output MeshBuilder
 * @param truss  input Truss
 * @param bColor Flag indicating whether to color the truss mesh should be calculated from hash of the block index.
 * @param width  The width of the truss lines or tubes. Default value is -1.
 */
void drawTruss_mesh( Mesh::Builder& mesh, const Truss& truss, bool bColor, float width=-1 ){
    //printf("=============\n");
    for( int i=0; i<truss.blocks.size(); i++ ){
        //Vec3f color;
        //if(bColor) Draw::color_of_hash(i+15454, mesh.penColor);
        Vec2i bj,bi = truss.blocks[i];
        if( i==(truss.blocks.size()-1) ){ bj = {truss.points.size(),truss.edges.size()}; }else{ bj = truss.blocks[i+1]; };
        //Draw3D::drawPoints( bj.a-bi.a, &truss.points[bi.a], 0.1 );
        if(width>0){ mesh.newSub( Mesh::TRIANGLES ); }else{ mesh.newSub( Mesh::LINES ); }
        for( int j=bi.b; j<bj.b; j++ ){
            //Vec3f a,b;
            const TrussEdge& ed = truss.edges[j];
            Vec3f p1 = (Vec3f)truss.points[ed.a];
            Vec3f p2 = (Vec3f)truss.points[ed.b];
            if(width>0){
                Vec3f up,lf; (p2-p1).getSomeOrtho( up, lf );
                mesh.addTube4( p1, p2, up, width, width );
            }else{
                mesh.addLine( p1, p2 );
            }
            //if(i==(truss.blocks.size()-1)){ printf( "%i %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, ed.a, ed.b, a.x, a.y, a.z,  b.x, b.y, b.z ); }
        }
    }
    //printf( "points.size() %i \n", truss.points.size() );
};

/**
 * Draws a plate into MeshBuilder.     
 * Deprecated, use drawPlateContour_mesh() instead.
 *
 * @param mesh    output MeshBuilder
 * @param o       input SpaceCrafting::Plate to draw
 * @param nodes   array of SpaceCrafting::Node objects on which the plate girder are attached
 * @param girders array of SpaceCrafting::Girder objects on which the plate is attached
 */
void drawPlate_mesh( Mesh::Builder& mesh, const Plate& o, const Node* nodes, const Girder* girders ){
    Vec3f p00=(Vec3f)nodes[ girders[o.g1].nodes.x ].pos;
    Vec3f p01=(Vec3f)nodes[ girders[o.g1].nodes.y ].pos;
    Vec3f p10=(Vec3f)nodes[ girders[o.g2].nodes.x ].pos;
    Vec3f p11=(Vec3f)nodes[ girders[o.g2].nodes.y ].pos;
    Vec3f d;
    d = p01-p00; p01=p00+d*o.g1span.x;  p00.add_mul( d,o.g1span.y);
    d = p11-p10; p11=p10+d*o.g2span.x;  p10.add_mul( d,o.g2span.y);
    mesh.addQuad( p00, p01, p10, p11 );
}

/**
 * Draws a plate into MeshBuilder.
 *
 * @param mesh    output MeshBuilder
 * @param o       input SpaceCrafting::Plate to draw
 * @param nodes   array of SpaceCrafting::Node objects on which the plate girder are attached
 * @param girders array of SpaceCrafting::Girder objects on which the plate is attached
 * @param filled  flag indicating whether the plate should be filled or just the contour should be drawn
 */
void drawPlateContour_mesh( Mesh::Builder& mesh,  const Plate& o, const Node* nodes, const Girder* girders, bool filled ){
    Quad3d qd;
    qd.l1.fromSubLine( nodes[ girders[o.g1].nodes.x ].pos, nodes[ girders[o.g1].nodes.y ].pos, o.g1span.x, o.g1span.y );
    qd.l2.fromSubLine( nodes[ girders[o.g2].nodes.x ].pos, nodes[ girders[o.g2].nodes.y ].pos, o.g2span.x, o.g2span.y );
    //printf( );
    //Draw3D::drawQuad(qd, filled);
    mesh.addQuad( (Vec3f)qd.p00, (Vec3f)qd.p01, (Vec3f)qd.p10, (Vec3f)qd.p11 );
    //Draw3D::drawQuad_bare();
}

/**
 * Draws SpaceCraft into MeshBuilder. Level-of-detail (LOD) can be specified. If LOD is 0, the mesh will be drawn as lines. If LOD is 1, the mesh will be drawn as triangles.
 *
 * @param spaceCraft input SpaceCraft object to be drawn.
 * @param mesh       output MeshBuilder object to store the generated mesh.
 * @param iLOD       level of detail for the mesh generation.
 * @param bText      Flag indicating whether to include text in the mesh. ( Not implemented yet )
 * @param bColor     Flag indicating whether to use schematic color for the mesh rendering.
 * @param color      The color to be used for the mesh.  (cuttently not used)
 * @return           The number of submeshes in the generated mesh.
 */
int drawSpaceCraft_Mesh( const SpaceCraft& spaceCraft, Mesh::Builder& mesh, int iLOD, bool bText, bool bColor, Vec3f color ){
    //printf( "##### START drawSpaceCraft_Mesh() \n" );
    const std::vector<Node>&   nodes   = spaceCraft.nodes;
    const std::vector<Girder>& girders = spaceCraft.girders;
    // --- Ropes
    float line_width = 0.25;   // TODO : this should be set better
    if( iLOD>0 ){
        if(bColor) mesh.penColor = (Vec3f){0.0,0.0,0.0};
        drawTruss_mesh( mesh, spaceCraft.truss, false, line_width );
    }
    if(!bColor) mesh.penColor = (Vec3f){0.2,0.2,0.2};
    for( const Rope& rp : spaceCraft.ropes ){
        Vec3f p1=(Vec3f)nodes[rp.nodes.x].pos;
        Vec3f p2=(Vec3f)nodes[rp.nodes.y].pos;
        if(iLOD==0){
            mesh.newSub( Mesh::LINES ); // toDo
            if(line_width>0){
                Vec3f up,lf; (p2-p1).getSomeOrtho( up, lf );
                mesh.addTube4( p1, p2, up, line_width, line_width );
            }else{
                mesh.addLine( p1, p2 );
            }
        }
        // ToDo : polygonized line (tube) with non-zero width ?
    };
    // --- Girders
    if(!bColor) mesh.penColor = (Vec3f){0.1,0.1,0.5};
    for( const Girder& gd : spaceCraft.girders ){
        Vec3f p0=(Vec3f)nodes[gd.nodes.x].pos;
        Vec3f p1=(Vec3f)nodes[gd.nodes.y].pos;
        if(iLOD==0){
            mesh.newSub( Mesh::LINES ); // toDo
            mesh.addLine( p0, p1 );
        }
    };
    for( const Gun& o : spaceCraft.guns ){
        Vec3f p0=(Vec3f)nodes[ girders[o.suppId].nodes.x ].pos;
        Vec3f p1=(Vec3f)nodes[ girders[o.suppId].nodes.y ].pos;
        Vec3f d;
        d = p1-p0; p1=p0+d*o.suppSpan.x;  p0.add_mul( d,o.suppSpan.y);
        if(iLOD==0){
            mesh.newSub( Mesh::LINES );
        }
        if(iLOD>0){
            mesh.penColor = (Vec3f){0.9,0.9,0.9};
            //mesh.newSub(GL_TRIANGLES);
            mesh.newSub(Mesh::TRIANGLE_STRIP);
            //Mesh::HarmonicTube( {200,8}, {0.0,0.0},{1.0,2*M_PI}, 2.0, 2.0, d.norm(), 0.0, 100.0*2*M_PI, 0.5, false );
            Mesh::HarmonicTube2Mesh      ( {200,8}, {0.0,0.0},{1.0,2*M_PI}, 2.0, 2.0, d.norm(), 0.0, 100.0*2*M_PI, 0.5, false, mesh );
            mesh.move( mesh.subVertRange(-1), p1  );
            //glPopMatrix();
        }
    };
    // --- Rings
    for( const Ring& o : spaceCraft.rings ){
        if(iLOD==0){
            mesh.newSub(Mesh::LINES);
            mesh.addCircleAxis( o.nseg, (Vec3f)o.pose.pos, (Vec3f)o.pose.rot.b, (Vec3f)o.pose.rot.c, o.R );
        };
    };
    // --- Radiators
    for( const Radiator& o : spaceCraft.radiators ){
        drawPlate_mesh( mesh, o, nodes.data(), girders.data());
    };
    // --- Shields
    for( const Shield& o : spaceCraft.shields ){
        drawPlate_mesh(mesh, o, nodes.data(), girders.data() );
    };
    // --- Tanks
    for( const Tank& o : spaceCraft.tanks ){
        //Draw3D::drawCapsula( (Vec3f)(o.pose.pos+o.pose.rot.c*(o.span.c*0.5)), Vec3f(o.pose.pos+o.pose.rot.c*(o.span.c*-0.5)), o.span.a, o.span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        mesh.newSub(Mesh::TRIANGLES);
        Mesh::drawCapsula( mesh, (Vec3f)(o.pose.pos+o.pose.rot.c*(o.span.c*0.5)), Vec3f(o.pose.pos+o.pose.rot.c*(o.span.c*-0.5)), o.span.a, o.span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
    };
    for( const Thruster& o : spaceCraft.thrusters ){
        float R = o.span.b;
        float L = o.span.c;
        if( iLOD>0 ){
            mesh.newSub( Mesh::TRIANGLE_STRIP );
            Parabola2Mesh_ExtrudedWire ( {8,16}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, 5.0  , mesh  );
        }else{
            mesh.newSub( Mesh::TRIANGLE_STRIP );
            Parabola2Mesh              ( {16,16}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, true, mesh );
        }
        mesh.applyMatrix( mesh.subVertRange(-1), (Mat3f)o.pose.rot );
        mesh.move       ( mesh.subVertRange(-1), (Vec3f)o.pose.pos );
    };
    /*
    for( const Rock& o : spaceCraft.rocks ){
        glPushMatrix();
        //printf( "rock.pose.rot \n"); o.pose.rot.print(); printf( "Rock pos(%f,%f,%f) span(%f,%f,%f) \n", o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,   o.span.x, o.span.y, o.span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o.pose.pos, o.pose.rot, o.span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_asteroide);
        glPopMatrix();
    }
    */
    /*
    for( const Balloon& o : spaceCraft.balloons ){
        glPushMatrix();
        printf( "ballon.pose.rot \n"); o.pose.rot.print(); printf( "ballon pos(%f,%f,%f) span(%f,%f,%f) \n", o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,   o.span.x, o.span.y, o.span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o.pose.pos, o.pose.rot, o.span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_geoSphere);
        glPopMatrix();
    }
    */
    // --- Guns

    //for( Rope& rp : ropes ){ drawRope(rp);
    //
    //};
    //glLineWidth(1);
    printf( "##### END drawSpaceCraft_Mesh() \n" );
    return 0;
};


} // namespace SpaceCrafting

#endif
