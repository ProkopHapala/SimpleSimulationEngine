
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

#include "Noise.h"
#include "SphereSampling.h"
#include "DrawSphereMap.h"

#include "SpaceCraft.h"
#include "EditSpaceCraft.h"

namespace SpaceCrafting{

int ogl_asteroide=0, ogl_geoSphere=0;

void drawAsteroide( int nsamp, int nCrater, float hscale, bool wire ){

    //const int nsamp = 32;
    //const int nsamp = 8;
    const int npix  = 10*nsamp*nsamp;
    //uint32_t pix   [npix];
    float*  heights = new float[npix];
    //Vec3f samplePs[npix];
    Vec3d*  craterPos = new Vec3d[nCrater];
    double* craterSz  = new double[nCrater];

    //srand(34646);
    for(int i=0; i<nCrater; i++){
        //craterPos[i].setHomogenousSphericalSample( randf(), randf() );
        craterPos[i].fromRandomSphereSample();
        //printf( "%i %f %f %f \n", craterPos[i].x, craterPos[i].y, craterPos[i].z);
        craterSz[i] = randf()+0.4;
    }

    Vec2d dab = (Vec2d){ 1.0/nsamp, 1.0/nsamp };
    for(int iface=0; iface<10; iface++ ){
        for(int ia=0; ia<nsamp; ia++ ){
            for(int ib=0; ib<nsamp; ib++ ){
                Vec3d p;
                SphereSampling::icosa2cartes( (Vec2i){nsamp,nsamp}, iface, ia*dab.a, ib*dab.b, p );
                int i = iface*nsamp*nsamp + ia*nsamp + ib;

                //heights[i] = getSphericalHarmonicRand( {5,4}, p );

                p.normalize();
                double h;
                h = Noise::getQuadruRand( p, 35456 )*2;

                h += Noise::getCraterHeight( p, nCrater, 1.0, craterPos, craterSz )*0.3;
                heights[i] = h; //+ sin( h*4 );

            }
        }
    };


    //shape=glGenLists(1);
	//glNewList( shape, GL_COMPILE );
        //printf( " Solids::nTetrahedron_tris %i \n", Solids::nTetrahedron_tris );
        //for(int i=0; i<nCrater; i++){ Draw3D::drawPointCross( craterPos[i], craterSz[i] );}
    glDisable(GL_LIGHTING);
    glShadeModel( GL_SMOOTH     );

    if(wire) SphereSampling::drawIcosaMapWire( (Vec2i){nsamp,nsamp}, heights, hscale );
    else     SphereSampling::drawIcosaMap    ( (Vec2i){nsamp,nsamp}, heights, hscale );
    //glPopMatrix();
	//glEndList();

    //float heights = float[npix];
    //Vec3f samplePs[npix];
    //Vec3d  craterPos = new Vec3d[nCrater];
    //double craterSz  = new double[nCrater];

	delete [] heights;
	delete [] craterPos;
	delete [] craterSz;

}

void drawTruss( const Truss& truss ){
    //printf("=============\n");
    for( int i=0; i<truss.blocks.size(); i++ ){
        Draw::color_of_hash(i+15454);
        Vec2i bj,bi = truss.blocks[i];
        if( i==(truss.blocks.size()-1) ){ bj = {truss.points.size(),truss.edges.size()}; }else{ bj = truss.blocks[i+1]; };
        //Draw3D::drawPoints( bj.a-bi.a, &truss.points[bi.a], 0.1 );
        glBegin( GL_LINES );
        for( int j=bi.b; j<bj.b; j++ ){
            Vec3f a,b;
            const TrussEdge& ed = truss.edges[j];
            convert( truss.points[ed.a], a ); glVertex3f( a.x, a.y, a.z );
            convert( truss.points[ed.b], b ); glVertex3f( b.x, b.y, b.z );
            //if(i==(truss.blocks.size()-1)){ printf( "%i %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, ed.a, ed.b, a.x, a.y, a.z,  b.x, b.y, b.z ); }
        }
        glEnd();
    }
    //printf( "points.size() %i \n", truss.points.size() );
};

void drawPlate( const Plate& o, const Node* nodes, const Girder* girders ){
    Vec3f p00=(Vec3f)nodes[ girders[o.g1].p0 ].pos;
    Vec3f p01=(Vec3f)nodes[ girders[o.g1].p1 ].pos;
    Vec3f p10=(Vec3f)nodes[ girders[o.g2].p0 ].pos;
    Vec3f p11=(Vec3f)nodes[ girders[o.g2].p1 ].pos;
    Vec3f d;
    d = p01-p00; p01=p00+d*o.g1span.x;  p00.add_mul( d,o.g1span.y);
    d = p11-p10; p11=p10+d*o.g2span.x;  p10.add_mul( d,o.g2span.y);
    //print(p00); print(p01); print(p10); print(p11); printf("\n");
    Vec3f nor; nor.set_cross( p10-p00, p01-p00 ); nor.normalize();
    //Draw3D::drawVecInPos(nor*100.0,p00);
    glBegin(GL_QUADS);
    glNormal3f(nor.x,nor.y,nor.z);
    glVertex3f(p00.x,p00.y,p00.z);
    glVertex3f(p01.x,p01.y,p01.z);
    glVertex3f(p11.x,p11.y,p11.z);
    glVertex3f(p10.x,p10.y,p10.z);
    glEnd();
    //sprintf( str, "Radiator_%03i", gd.id );
    //Vec3f d = p1-p0; d.normalize();
    //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
    //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
}

void drawSpaceCraft( const SpaceCraft& spaceCraft, int iLOD ){
    //nodes,ropes;girders;thrustes;guns;radiators;shields;tanks;pipes;
    //for( None nd : nodes ){};
    const std::vector<Node>& nodes   = spaceCraft.nodes;
    const std::vector<Girder>& girders = spaceCraft.girders;
    // --- Ropes
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    if( iLOD>0 ){
        glLineWidth(1.0);
        glColor3f(0.0,0.0,0.0);
        drawTruss( spaceCraft.truss );
    }

    glLineWidth(0.5);
    glColor3f(0.0,0.0,0.0);
    for( const Rope& rp : spaceCraft.ropes ){
        Vec3f p0=(Vec3f)nodes[rp.p0].pos;
        Vec3f p1=(Vec3f)nodes[rp.p1].pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
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
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
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
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
        if(iLOD>0){
            glPushMatrix();
            glTranslatef( p1.x, p1.y, p1.z );
            glEnable( GL_DEPTH_TEST );
            glEnable( GL_LIGHTING   );
            glColor3f(0.9,0.9,0.9);
            Draw3D::drawUV_HarmonicTube( {200,8}, {0.0,0.0},{1.0,2*M_PI}, 2.0, 2.0, d.norm(), 0.0, 100.0*2*M_PI, 0.5, false );
            glPopMatrix();
        }
    };
    // --- Rings
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    glColor3f(0.1,0.1,0.5);
    for( const Ring& o : spaceCraft.rings ){
        //Vec3f p0=(Vec3f)nodes[o.p0].pos;
        //Vec3f p1=(Vec3f)nodes[o.p1].pos;
        //Vec3f uax = (p1-p0); uax.normalize();
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        //Draw3D::drawCircleAxis( o.nseg, (p0+p1)*0.5f, (Vec3f)o.up, uax, o.R );
        if(iLOD==0) Draw3D::drawCircleAxis( o.nseg, o.pose.pos, o.pose.rot.b, o.pose.rot.c, o.R );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );
        //sprintf( str, "Girder_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.p0].pos+nodes[rp.p1].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Radiators
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glColor3f(0.3,0.1,0.1); // ToDo by temperature
    for( const Radiator& o : spaceCraft.radiators ){
        drawPlate(o, nodes.data(), girders.data());
    };
    // --- Shields
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glColor3f(0.9,0.9,0.9);
    for( const Shield& o : spaceCraft.shields ){
        drawPlate(o, nodes.data(), girders.data() );
    };
    // --- Tanks
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glColor3f(0.5,0.5,0.5);
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
    glColor3f(0.5,0.5,0.5);
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
        //Draw3D::drawUV_Parabola( {16,16}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, true );
        if( iLOD>0 ){ Draw3D::drawUV_Parabola_ExtrudedWire( {8,16}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, 5.0 ); }
        else        { Draw3D::drawUV_Parabola( {16,16}, {0.0,0.0}, {1.0,M_PI*2}, R, -L, 0.5, true );  }
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
    for( const Rock& o : spaceCraft.rocks ){
        glPushMatrix();
        //printf( "rock.pose.rot \n"); o.pose.rot.print(); printf( "Rock pos(%f,%f,%f) span(%f,%f,%f) \n", o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,   o.span.x, o.span.y, o.span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o.pose.pos, o.pose.rot, o.span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_asteroide);
        glPopMatrix();
    }

    for( const Balloon& o : spaceCraft.balloons ){
        glPushMatrix();
        printf( "ballon.pose.rot \n"); o.pose.rot.print(); printf( "ballon pos(%f,%f,%f) span(%f,%f,%f) \n", o.pose.pos.x, o.pose.pos.y, o.pose.pos.z,   o.span.x, o.span.y, o.span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o.pose.pos, o.pose.rot, o.span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_geoSphere);
        glPopMatrix();
    }

    // --- Guns

    //for( Rope& rp : ropes ){ drawRope(rp);
    //
    //};
    glLineWidth(1);
};





} // namespace SpaceCrafting

#endif
