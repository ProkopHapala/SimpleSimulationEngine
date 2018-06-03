
#ifndef  SpaceCraftDraw_h
#define  SpaceCraftDraw_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "Draw3D_Surf.h"
#include "SDL_utils.h"

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "SpaceCraft.h"
#include "SpaceCraft.h"

#include "SpaceCraft.h"
#include "EditSpaceCraft.h"

namespace SpaceCrafting{

void drawTruss( Truss& truss ){
    //printf("=============\n");
    for( int i=0; i<truss.blocks.size(); i++ ){
        Draw::color_of_hash(i+15454);
        Vec2i bj,bi = truss.blocks[i];
        if( i==(truss.blocks.size()-1) ){ bj = {truss.points.size(),truss.edges.size()}; }else{ bj = truss.blocks[i+1]; };
        //Draw3D::drawPoints( bj.a-bi.a, &truss.points[bi.a], 0.1 );
        glBegin( GL_LINES );
        for( int j=bi.b; j<bj.b; j++ ){
            Vec3f a,b;
            TrussEdge& ed = truss.edges[j];
            convert( truss.points[ed.a], a ); glVertex3f( a.x, a.y, a.z );
            convert( truss.points[ed.b], b ); glVertex3f( b.x, b.y, b.z );

            //if(i==(truss.blocks.size()-1)){ printf( "%i %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, ed.a, ed.b, a.x, a.y, a.z,  b.x, b.y, b.z ); }
        }
        glEnd();
    }
    //printf( "points.size() %i \n", truss.points.size() );
};

void drawSpaceCraft( const SpaceCraft& spaceCraft ){
    //nodes,ropes;girders;thrustes;guns;radiators;shields;tanks;pipes;
    //for( None nd : nodes ){};
    const std::vector<Node>& nodes   = spaceCraft.nodes;
    const std::vector<Girder>& girders = spaceCraft.girders;
    // --- Ropes
    glDisable(GL_CULL_FACE);
    glLineWidth(0.5);
    glColor3f(0.0,0.0,0.0);
    for( const Rope& rp : spaceCraft.ropes ){
        Vec3f p0=(Vec3f)nodes[rp.p0].pos;
        Vec3f p1=(Vec3f)nodes[rp.p1].pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        Draw3D::drawLine( p0,p1 );
        sprintf( str, "rope_%03i", rp.id );
        Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, (Vec3f){1.0,0.0,0.0}, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Girders
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    glColor3f(0.1,0.1,0.5);
    for( const Girder& gd : spaceCraft.girders ){
        Vec3f p0=(Vec3f)nodes[gd.p0].pos;
        Vec3f p1=(Vec3f)nodes[gd.p1].pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        Draw3D::drawLine( p0,p1 );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );

        sprintf( str, "Girder_%03i", gd.id );
        Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    glLineWidth(5);
    glColor3f(0.6,0.1,0.1);
    for( const Gun& o : spaceCraft.guns ){
        Vec3f p0=(Vec3f)nodes[ girders[o.suppId].p0 ].pos;
        Vec3f p1=(Vec3f)nodes[ girders[o.suppId].p1 ].pos;
        Vec3f d;
        d = p1-p0; p1=p0+d*o.suppSpan.x;  p0.add_mul( d,o.suppSpan.y);
        Draw3D::drawLine( p0,p1 );
    };
    // --- Rings
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    glColor3f(0.1,0.1,0.5);
    for( const Ring& o : spaceCraft.rings ){
        Vec3f p0=(Vec3f)nodes[o.p0].pos;
        Vec3f p1=(Vec3f)nodes[o.p1].pos;
        Vec3f uax = (p1-p0); uax.normalize();
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        Draw3D::drawCircleAxis( o.nseg, (p0+p1)*0.5f, (Vec3f)o.up, uax, o.R );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );

        //sprintf( str, "Girder_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Radiators
    glLineWidth(3);
    printf("");
    glColor3f(0.5,0.5,0.5);
    //glEnable( GL_BLEND );
    glEnable(GL_DEPTH_TEST);
    for( const Radiator& rd : spaceCraft.radiators ){
        Vec3f p00=(Vec3f)nodes[ girders[rd.g1].p0 ].pos;
        Vec3f p01=(Vec3f)nodes[ girders[rd.g1].p1 ].pos;
        Vec3f p10=(Vec3f)nodes[ girders[rd.g2].p0 ].pos;
        Vec3f p11=(Vec3f)nodes[ girders[rd.g2].p1 ].pos;
        Vec3f d;
        d = p01-p00; p01=p00+d*rd.g1span.x;  p00.add_mul( d,rd.g1span.y);
        d = p11-p10; p11=p10+d*rd.g2span.x;  p10.add_mul( d,rd.g2span.y);
        //print(p00); print(p01); print(p10); print(p11); printf("\n");
        glBegin(GL_QUADS);
        glVertex3f(p00.x,p00.y,p00.z);
        glVertex3f(p01.x,p01.y,p01.z);
        glVertex3f(p11.x,p11.y,p11.z);
        glVertex3f(p10.x,p10.y,p10.z);
        glEnd();
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Tanks
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    for( const Tank& o : spaceCraft.tanks ){
        Draw3D::drawCapsula( (Vec3f)(o.pose.pos+o.pose.rot.c*(o.span.c*0.5)), Vec3f(o.pose.pos+o.pose.rot.c*(o.span.c*-0.5)), o.span.a, o.span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };

    // --- Thrusters
    glLineWidth(1);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    for( const Thruster& o : spaceCraft.thrusters ){

        //Draw3D::drawMatInPos( (Mat3f)o.pose.rot*10.0, (Vec3f)o.pose.pos );
        glPushMatrix();
        //glTranslatef( o.pose.pos.x, o.pose.pos.y, o.pose.pos.z );
        Draw3D::rigidTransform( (Vec3f)o.pose.pos, (Mat3f)o.pose.rot, (Vec3f){1.0,1.0,1.0}, false);
        float R = o.span.b;
        float L = o.span.c;
        printf( "Thruster R %f L %f", R, L );
    //  drawUV_Cone( Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire ){
        //Draw3D::drawUV_Cone( {10,10}, {0.0,0.0}, {1.0,M_PI*2}, o.span.a, o.span.b, o.span.c, 0.5, false );
        //Draw3D::drawUV_Cone( {10,10}, {0.0,0.0}, {1.0,M_PI*2}, o.span.a, o.span.b, o.span.c, 0.5, true );
        //Draw3D::drawUV_Teardrop( {10,10}, {0.0,0.0}, {1.0,M_PI*2}, o.span.a, o.span.b*0.1, o.span.c, 0.5, true );
        Draw3D::drawUV_Parabola( {10,10}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, true );
        //glColor3f(1.0,0.0,0.0);
        //Draw3D::drawUV_Hyperbola( {20,20}, {0.0,0.0}, {1.0,M_PI*2},  o.span.a, o.span.b, o.span.c, 0.5, true );
        //Draw3D::drawUV_Hyperbola( {20,20}, {0.0,0.0}, {1.0,M_PI*2},   5.0, 8.0, 5.0, 0.5, true );
        //glColor3f(0.0,0.0,1.0);
        //Draw3D::drawUV_Hyperbola( {20,20}, {0.0,0.0}, {1.0,M_PI*2}, -o.span.a, o.span.b, o.span.c, 0.5, true );
        //Draw3D::drawUV_Hyperbola( {20,20}, {0.0,0.0}, {1.0,M_PI*2},  -5.0, 8.0, 5.0, 0.5, true );

        //Draw3D::drawShere( (Vec3f)(o.pose.pos+o.pose.rot.c*(o.span.c*0.5)), Vec3f(o.pose.pos+o.pose.rot.c*(o.span.c*-0.5)), o.span.a, o.span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        glPopMatrix();
    };

    // --- Guns

    //for( Rope& rp : ropes ){ drawRope(rp);
    //
    //};
    glLineWidth(1);
};





} // namespace SpaceCrafting

#endif
