
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
#include "MeshBuilder.h"

namespace SpaceCrafting{

int ogl_asteroide=0, ogl_geoSphere=0;

void renderTruss(int nb, int2* bonds, Quat4f* ps, float* strain=0, float sc=1.0 ){
    glBegin(GL_LINES);
    for(int i=0; i<nb; i++ ){
        //printf( "renderTruss()[%i] \n", i );
        int2 b =  bonds[i];
        if(strain){
            float f=strain[i];
            //if( fabs(f)>0.0001 )printf( "Edge[%i] strain=%g \n", i, f );
            f*=sc;
            if(f>0){ 
                //Draw3D::color(Vec3f{f*0,0,0}); 
                glColor3f(f,f*0.2,0.);
            }else{ 
                //Draw3D::color(Vec3f{0,f*100.0f,f});
                glColor3f(0.,f*-0.2,-f);  
                //glColor3f(0,1.,0.);
                //glColor3f(0,-f*10000000.0,0.);
                //printf("neg\n"); 
            };
            //Draw3D::color(Vec3f{0.,1.0,0.});
        } 
        Draw3D::vertex( ps[b.x].f );
        Draw3D::vertex( ps[b.y].f );
        //printf( "renderTruss[%i](%i,%i) p(%g,%g,%g) p(%g,%g,%g)\n", i, b.x, b.y,  ps[b.x].f.x,ps[b.x].f.y,ps[b.x].f.z,   ps[b.y].f.x,ps[b.y].f.y,ps[b.y].f.z );
    }
    glEnd();
}

void renderPoinSizes(int n, Quat4f* ps, float sc=1.0 ){
    glColor3f(0.0,0.0,1.0);
    glPointSize( 10.0 );
    glBegin(GL_POINTS);
    for(int i=0; i<n; i++ ){
        //printf( "point[%i] m=%g \n", i, ps[i].e );
        float c = ps[i].e*sc;
        glColor3f(c,c,c);
        //glPointSize( ps[i].e*sc );
        Draw3D::vertex( ps[i].f );
    }
    glEnd();
    //exit(0);
}

void renderPointForces(int n, Quat4f* ps, Quat4f* fs, float sc=1.0 ){
    glColor3f(1.0,0.0,1.0);
    //glPointSize( 10.0 );
    glBegin(GL_LINES);
    for(int i=0; i<n; i++ ){
        //printf( "point[%i] m=%g \n", i, ps[i].e );
        //float c = ps[i].e*sc;
        //glColor3f(c,c,c);
        //glPointSize( ps[i].e*sc );
        Draw3D::vertex( ps[i].f );
        Draw3D::vertex( ps[i].f+fs[i].f*sc );
    }
    glEnd();
    //exit(0);
}

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


void drawPlate( const Plate& o, const Girder*const* girders ){
    Vec3f p00=(Vec3f)girders[o.g1]->nodes.x->pos;
    Vec3f p01=(Vec3f)girders[o.g1]->nodes.y->pos;
    Vec3f p10=(Vec3f)girders[o.g2]->nodes.x->pos;
    Vec3f p11=(Vec3f)girders[o.g2]->nodes.y->pos;
    Vec3f d;
    d = p01-p00; p01=p00+d*o.g1span.x;  p00.add_mul( d,o.g1span.y);
    d = p11-p10; p11=p10+d*o.g2span.x;  p10.add_mul( d,o.g2span.y);
    //print(p00); print(p01); print(p10); print(p11); printf("\n");
    Vec3f nor; nor.set_cross( p10-p00, p01-p00 ); nor.normalize();
    //Draw3D::drawVecInPos(nor*100.0,p00);
    //mesh.addQuad( p00, p01, p10, p11 );
    glBegin(GL_QUADS);
    glNormal3f(nor.x,nor.y,nor.z);
    glVertex3f(p00.x,p00.y,p00.z);
    glVertex3f(p01.x,p01.y,p01.z);
    glVertex3f(p11.x,p11.y,p11.z);
    glVertex3f(p10.x,p10.y,p10.z);
    glEnd();
}

void drawPlateContour( const Plate& o, const Girder*const* girders, bool filled ){
    Quad3d qd;
    qd.l1.fromSubLine( girders[o.g1]->nodes.x ->pos, girders[o.g1]->nodes.y->pos, o.g1span.x, o.g1span.y );
    qd.l2.fromSubLine( girders[o.g2]->nodes.x ->pos, girders[o.g2]->nodes.y->pos, o.g2span.x, o.g2span.y );
    //printf( );
    Draw3D::drawQuad(qd, filled);
    //mesh.addQuad( qd.p00, qd.p01, qd.p10, qd.p11 );
    //Draw3D::drawQuad_bare();
}

int drawSliderPath( const Slider* o, const Quat4f* ps, float sz=10 ){
    Draw3D::drawPointCross( ps[o->along.x].f, sz );
    Draw3D::drawLineStrip( o->path.n, o->path.ps, ps, o->path.closed );
    return 0;
}

/*
int drawSpaceCraft_sliderPaths( const SpaceCraft& craft, const Quat4f* ps, float sz=10 ){
    //Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    int n = craft.sliders.size();
    for(int i=0; i<n; i++){
        const Slider* o = craft.sliders[i];
        Draw3D::drawPointCross( ps[o->ifix].f, sz );
        Draw3D::drawLineStrip( o->path.n, o->path.ps, ps, o->path.closed );
    }
    return n;
}
*/

void drawNode( const Node& o, float sz=10 ){
    //Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    Draw3D::drawPointCross( o.pos, sz );
}

void drawNode( const Node& o, const Quat4f* ps, float sz=10.f ){
    //Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    if(o.boundTo){
        Draw3D::drawPointCross( ps[o.along.x+o.boundTo->pointRange.x].f, sz );
    }else{
        Draw3D::drawPointCross( ps[o.id].f, sz );
    }
}

void drawSpaceCraft_nodes( const SpaceCraft& craft, const Quat4f* ps=0, float sz=10, int i0=0, int i1=-1 ){
    //Draw3D::color(Vec3f{1.0f,0.0f,1.0f});
    if(i1<0) i1 = craft.nodes.size() + i1;
    if( ps ){ for(int i=i0; i<=i1; i++){ drawNode( *craft.nodes[i], ps, sz ); } }
    else    { for(int i=i0; i<=i1; i++){ drawNode( *craft.nodes[i],     sz ); } }
}

int drawPointRange( int n, Vec2i prange, int byN, int off, Vec2d span, const Quat4f* ps ){
    //printf( "drawPointRange(n=%i prange=%i byN=%i off=%i span(%g,%g))\n", n, prange.x, prange.y, byN, off, span.x, span.y );
    int i0 = prange.x;
    int np = (prange.y - prange.x)/byN;
    //double d = (span.y-span.x)/n;
    double step = 1./(n-1);
    glBegin(GL_POINTS);
    for(int i=0;i<n; i++){
        double c = i*step;
        int ip = prange.x + off + byN*(int)( ( span.x*(1-c) + span.y*c )*np  + 0.5 );
        //int ip = i0 + off + byN*i;
        //printf( "plateOnGriders()[%i] (%i/%i) (%i/%i)\n", i, i1,n1, i2,n2 );
        Draw3D::vertex( ps[ip].f );      
    }
    /*
    for(int i=0;i<np; i++){
        int ip = i0 + i*byN + off;
        //printf( "plateOnGriders()[%i] (%i/%i) (%i/%i)\n", i, i1,n1, i2,n2 );
        if( (i>(span.x*np))&&(i<(span.y*np)) ){
            Draw3D::vertex( ps[ip].f );
        }      
    }
    */
    //return ibloc;
    glEnd();
    return 0;
}

void drawPicked( const SpaceCraft& craft, int ipick ){
    switch( (ComponetKind)craft.pickedTyp){
        case ComponetKind::Radiator :{
            //printf("drawPicked radiator[%i] \n", ipick );
            drawPlateContour( *craft.radiators[ipick], &craft.girders[0], false );
            }break;
        case ComponetKind::Shield:{
            //printf("drawPicked shield[%i] \n", ipick );
            drawPlateContour( *craft.shields  [ipick], &craft.girders[0], false );
            }break;
        case ComponetKind::Girder:{
            const Girder& o = *craft.girders[ipick];
            Draw3D::drawLine(o.nodes.x->pos,o.nodes.y->pos);
            } break;
        case ComponetKind::Rope:{
            const Rope& o = *craft.ropes[ipick];
            Draw3D::drawLine(o.nodes.x->pos,o.nodes.y->pos);
            } break;
    }
}


void drawSpaceCraft( const SpaceCraft& craft, int iLOD, bool bText, bool bColor ){
    // --- Ropes
    glLineWidth(0.5);
    if(!bColor)glColor3f(0.2,0.2,0.2);
    for( const Rope* o : craft.ropes ){
        Vec3f p0=(Vec3f)o->nodes.x->pos;
        Vec3f p1=(Vec3f)o->nodes.y->pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        //if(iLOD==0) 
        Draw3D::drawLine( p0,p1 );
        // if(bText){
        //     sprintf( str, "rope_%03i", rp.id );
        //     Vec3f d = p1-p0; d.normalize();
        //     //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //     Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        //     //Draw3D::drawText3D( str, (p0+p1)*0.5, (Vec3f){1.0,0.0,0.0}, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        // }
    }
    // --- Ropes
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    //if(!bColor)
    glColor3f(0.4,0.4,0.4);
    // --- Girders
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.1,0.1,0.5);
    for( const Girder* o : craft.girders ){
        Vec3f p0=(Vec3f)o->nodes.x->pos;
        Vec3f p1=(Vec3f)o->nodes.y->pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        if(iLOD==0){
            Draw3D::drawLine( p0,p1 );
        }else{
            Draw3D::drawCylinderStrip( 4, 1.0,1.0, p0, p1 );
        } 
        /*
        if(bText){
            sprintf( str, "Girder_%03i", gd.id );
            Vec3f d = p1-p0; d.normalize();
            //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
            Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        }
        */
    };
    glLineWidth(5);
    if(bColor)glColor3f(0.6,0.1,0.1);
    for( const Gun* o : craft.guns ){
        Vec3f p0=(Vec3f)craft.girders[o->suppId]->nodes.x->pos;
        Vec3f p1=(Vec3f)craft.girders[o->suppId]->nodes.y->pos;
        Vec3f d;
        d = p1-p0; p1=p0+d*o->suppSpan.x;  p0.add_mul( d,o->suppSpan.y);
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
        if(iLOD>0){
            glPushMatrix();
            glTranslatef( p1.x, p1.y, p1.z );
            glEnable( GL_DEPTH_TEST );
            glEnable( GL_LIGHTING   );
            glDisable(GL_CULL_FACE);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
            glColor3f(0.9,0.9,0.9);
            Draw3D::drawUV_HarmonicTube( {200,8}, {0.0,0.0},{1.0,2*M_PI}, 2.0, 2.0, d.norm(), 0.0, 100.0*2*M_PI, 0.5, false );
            glPopMatrix();
        }
    };
    // --- Rings
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.1,0.1,0.5);
    for( const Ring* o : craft.rings ){
        //Vec3f p0=(Vec3f)o.nodes.x->pos;
        //Vec3f p1=(Vec3f)o.nodes.y->pos;
        //Vec3f uax = (p1-p0); uax.normalize();
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        //Draw3D::drawCircleAxis( o.nseg, (p0+p1)*0.5f, (Vec3f)o.up, uax, o.R );
        if(iLOD==0) Draw3D::drawCircleAxis( o->nseg, o->pose.pos, o->pose.rot.b, o->pose.rot.c, o->R );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );
        //sprintf( str, "Girder_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Radiators
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    if(bColor)glColor3f(0.3,0.1,0.1); // ToDo by temperature
    for( const Radiator* o : craft.radiators ){
        drawPlate(*o, craft.girders.data() );
    };
    // --- Shields
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.2,0.2,0.2);
    for( const Shield* o : craft.shields ){
        drawPlate(*o, craft.girders.data() );
    };
    
    // --- Tanks
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.5,0.5,0.5);
    for( const Tank* o : craft.tanks ){
        Draw3D::drawCapsula( (Vec3f)(o->pose.pos+o->pose.rot.c*(o->span.c*0.5)), Vec3f(o->pose.pos+o->pose.rot.c*(o->span.c*-0.5)), o->span.a, o->span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Thrusters
    glLineWidth(1);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.5,0.5,0.5);
    for( const Thruster* o : craft.thrusters ){
        //Draw3D::drawMatInPos( (Mat3f)o.pose.rot*10.0, (Vec3f)o.pose.pos );
        glPushMatrix();
        //glTranslatef( o.pose.pos.x, o.pose.pos.y, o.pose.pos.z );
        Draw3D::rigidTransform( (Vec3f)o->pose.pos, (Mat3f)o->pose.rot, (Vec3f){1.0,1.0,1.0}, false);
        float R = o->span.b;
        float L = o->span.c;
        //printf( "Thruster R %f L %f", R, L );
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

        //Draw3D::drawShere( (Vec3f)(o->pose.pos+o->pose.rot.c*(o->span.c*0.5)), Vec3f(o->pose.pos+o->pose.rot.c*(o->span.c*-0.5)), o->span.a, o->span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        glPopMatrix();
    };
    for( const Rock* o : craft.rocks ){
        glPushMatrix();
        //printf( "rock.pose.rot \n"); o->pose.rot.print(); printf( "Rock pos(%f,%f,%f) span(%f,%f,%f) \n", o->pose.pos.x, o->pose.pos.y, o->pose.pos.z,   o->span.x, o->span.y, o->span.z ); // print(o->pose.rot);
        Draw3D::rigidTransform( o->pose.pos, o->pose.rot, o->span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_asteroide);
        glPopMatrix();
    }
    for( const Balloon* o : craft.balloons ){
        glPushMatrix();
        printf( "ballon.pose.rot \n"); o->pose.rot.print(); printf( "ballon pos(%f,%f,%f) span(%f,%f,%f) \n", o->pose.pos.x, o->pose.pos.y, o->pose.pos.z,   o->span.x, o->span.y, o->span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o->pose.pos, o->pose.rot, o->span );
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

#ifdef Truss_h

void drawTrussDirect( Truss& truss ){
    glBegin(GL_LINES);
    for(int i=0; i<truss.edges.size(); i++){
        TrussEdge& ed = truss.edges[i];
        Draw::color_of_hash( ed.type*1545 + 456 );
        Draw3D::vertex(truss.points[ed.a]);
        Draw3D::vertex(truss.points[ed.b]);
    }
    glEnd();
}

void drawTruss( const Truss& truss, bool bColor ){
    //printf("=============\n");
    for( int i=0; i<truss.blocks.size(); i++ ){
        Vec3f color;
        if(bColor) Draw::color_of_hash(i+15454);

        Vec2i bj,bi = truss.blocks[i];
        if( i==(truss.blocks.size()-1) ){ bj = {truss.points.size(),truss.edges.size()}; }else{ bj = truss.blocks[i+1]; };
        //Draw3D::drawPoints( bj.a-bi.a, &truss.points[bi.a], 0.1 );
        //mesh.newSub( GL_LINES );
        glBegin(GL_LINES);
        for( int j=bi.b; j<bj.b; j++ ){
            Vec3f a,b;
            const TrussEdge& ed = truss.edges[j];
            convert( truss.points[ed.a], a ); glVertex3f( a.x, a.y, a.z );
            convert( truss.points[ed.b], b ); glVertex3f( b.x, b.y, b.z );
            //mesh.addLine( (Vec3f)truss.points[ed.a], (Vec3f)truss.points[ed.b], color );
            //if(i==(truss.blocks.size()-1)){ printf( "%i %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", i, j, ed.a, ed.b, a.x, a.y, a.z,  b.x, b.y, b.z ); }
        }
        glEnd();
    }
    //printf( "points.size() %i \n", truss.points.size() );
};

void drawSpaceCraft( const SpaceCraft& craft, const Truss& truss, int iLOD, bool bText, bool bColor ){
    //nodes,ropes;girders;thrustes;guns;radiators;shields;tanks;pipes;
    //for( None nd : nodes ){};
    //const std::vector<Node*>&   nodes  = craft.nodes;
    //const std::vector<Girder>& girders = craft.girders;

    // --- Ropes
    glLineWidth(0.5);
    if(!bColor)glColor3f(0.2,0.2,0.2);
    for( const Rope* o : craft.ropes ){
        Vec3f p0=(Vec3f)o->nodes.x->pos;
        Vec3f p1=(Vec3f)o->nodes.y->pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
        // if(bText){
        //     sprintf( str, "rope_%03i", rp.id );
        //     Vec3f d = p1-p0; d.normalize();
        //     //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //     Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        //     //Draw3D::drawText3D( str, (p0+p1)*0.5, (Vec3f){1.0,0.0,0.0}, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        // }
    }
    // --- Ropes
    glDisable(GL_CULL_FACE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    //if(!bColor)
    glColor3f(0.4,0.4,0.4);
    if( iLOD>0 ){
        glLineWidth(1.0);
        //if(bColor)glColor3f(0.0,0.0,0.0);
        //drawTruss( craft.truss, bColor );
        drawTruss( truss, false );
    }
    // --- Girders
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.1,0.1,0.5);
    for( const Girder* o : craft.girders ){
        Vec3f p0=(Vec3f)o->nodes.x->pos;
        Vec3f p1=(Vec3f)o->nodes.y->pos;
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );
        /*
        if(bText){
            sprintf( str, "Girder_%03i", gd.id );
            Vec3f d = p1-p0; d.normalize();
            //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
            Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        }
        */
    };
    glLineWidth(5);
    if(bColor)glColor3f(0.6,0.1,0.1);
    for( const Gun* o : craft.guns ){
        Vec3f p0=(Vec3f)craft.girders[o->suppId]->nodes.x->pos;
        Vec3f p1=(Vec3f)craft.girders[o->suppId]->nodes.y->pos;
        Vec3f d;
        d = p1-p0; p1=p0+d*o->suppSpan.x;  p0.add_mul( d,o->suppSpan.y);
        if(iLOD==0) Draw3D::drawLine( p0,p1 );
        if(iLOD>0){
            glPushMatrix();
            glTranslatef( p1.x, p1.y, p1.z );
            glEnable( GL_DEPTH_TEST );
            glEnable( GL_LIGHTING   );
            glDisable(GL_CULL_FACE);
            glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
            glColor3f(0.9,0.9,0.9);
            Draw3D::drawUV_HarmonicTube( {200,8}, {0.0,0.0},{1.0,2*M_PI}, 2.0, 2.0, d.norm(), 0.0, 100.0*2*M_PI, 0.5, false );
            glPopMatrix();
        }
    };
    // --- Rings
    glLineWidth(3);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.1,0.1,0.5);
    for( const Ring* o : craft.rings ){
        //Vec3f p0=(Vec3f)o.nodes.x->pos;
        //Vec3f p1=(Vec3f)o.nodes.y->pos;
        //Vec3f uax = (p1-p0); uax.normalize();
        glEnable( GL_BLEND );
        glEnable(GL_DEPTH_TEST);
        //Draw3D::drawCircleAxis( o.nseg, (p0+p1)*0.5f, (Vec3f)o.up, uax, o.R );
        if(iLOD==0) Draw3D::drawCircleAxis( o->nseg, o->pose.pos, o->pose.rot.b, o->pose.rot.c, o->R );
        //Draw3D::drawCylinderStrip( 6, 1.0,1.0, p0, p1 );
        //sprintf( str, "Girder_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Radiators
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    if(bColor)glColor3f(0.3,0.1,0.1); // ToDo by temperature
    for( const Radiator* o : craft.radiators ){
        drawPlate(*o, craft.girders.data() );
    };
    // --- Shields
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.2,0.2,0.2);
    for( const Shield* o : craft.shields ){
        drawPlate(*o, craft.girders.data() );
    };
    
    // --- Tanks
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.5,0.5,0.5);
    for( const Tank* o : craft.tanks ){
        Draw3D::drawCapsula( (Vec3f)(o->pose.pos+o->pose.rot.c*(o->span.c*0.5)), Vec3f(o->pose.pos+o->pose.rot.c*(o->span.c*-0.5)), o->span.a, o->span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
    };
    // --- Thrusters
    glLineWidth(1);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glDisable(GL_CULL_FACE);
    //glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    if(bColor)glColor3f(0.5,0.5,0.5);
    for( const Thruster* o : craft.thrusters ){
        //Draw3D::drawMatInPos( (Mat3f)o.pose.rot*10.0, (Vec3f)o.pose.pos );
        glPushMatrix();
        //glTranslatef( o.pose.pos.x, o.pose.pos.y, o.pose.pos.z );
        Draw3D::rigidTransform( (Vec3f)o->pose.pos, (Mat3f)o->pose.rot, (Vec3f){1.0,1.0,1.0}, false);
        float R = o->span.b;
        float L = o->span.c;
        //printf( "Thruster R %f L %f", R, L );
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

        //Draw3D::drawShere( (Vec3f)(o->pose.pos+o->pose.rot.c*(o->span.c*0.5)), Vec3f(o->pose.pos+o->pose.rot.c*(o->span.c*-0.5)), o->span.a, o->span.b, M_PI*0.5, M_PI*0.5, M_PI*0.1, 16, true );
        //sprintf( str, "Radiator_%03i", gd.id );
        //Vec3f d = p1-p0; d.normalize();
        //Draw3D::drawText( str, nodes[rp.nodes.x].pos+nodes[rp.nodes.y].pos, fontTex, 10.1, 0 );
        //Draw3D::drawText3D( str, (p0+p1)*0.5, d, (Vec3f){0.0,1.0,0.0},  fontTex, 3.0, 0 );
        glPopMatrix();
    };
    for( const Rock* o : craft.rocks ){
        glPushMatrix();
        //printf( "rock.pose.rot \n"); o->pose.rot.print(); printf( "Rock pos(%f,%f,%f) span(%f,%f,%f) \n", o->pose.pos.x, o->pose.pos.y, o->pose.pos.z,   o->span.x, o->span.y, o->span.z ); // print(o->pose.rot);
        Draw3D::rigidTransform( o->pose.pos, o->pose.rot, o->span );
        //glScalef(100.0,100.0,100.0);
        glCallList(ogl_asteroide);
        glPopMatrix();
    }
    for( const Balloon* o : craft.balloons ){
        glPushMatrix();
        printf( "ballon.pose.rot \n"); o->pose.rot.print(); printf( "ballon pos(%f,%f,%f) span(%f,%f,%f) \n", o->pose.pos.x, o->pose.pos.y, o->pose.pos.z,   o->span.x, o->span.y, o->span.z ); // print(o.pose.rot);
        Draw3D::rigidTransform( o->pose.pos, o->pose.rot, o->span );
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

#endif // Truss_h

} // namespace SpaceCrafting

#endif
