
#ifndef  SpaceDraw_h
#define  SpaceDraw_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "SpaceBodies.h"
#include "RublePile.h"

#include "DrawSphereMap.h"


namespace SpaceDraw{

Vec3d ref_pos = Vec3dZero;
Vec3d ref_vel = Vec3dZero;
double zoom   = 1.0;

float scF = 1.0; float scSz = 0.2;
double view_scale = 1/1e+9;
Vec3d  view_shift = (Vec3d){0.0,0.0,0.0};

int iTrjMin=0,iTrjMax=0;
bool bRefTrj = false;
int trjStep = 10;
SpaceBody* referenceBody = 0;

int trj_n = 0;

void asPoints( int n, const SpaceBody* bodies, double epoch ){
    glBegin(GL_POINTS);
    for(int i=0; i<n; i++){
        const SpaceBody& body = bodies[i];
        Vec3d p = body.pointAtEpoch(epoch);
        if(isnan(p.x)) continue;
        glVertex3f((p.x-ref_pos.x)*zoom,(p.y-ref_pos.x)*zoom,(p.z-ref_pos.x)*zoom );
    }
    glEnd();
}

void asCrosses( int n, const SpaceBody* bodies, double epoch, double sz ){
    for(int i=0; i<n; i++){
        const SpaceBody& body = bodies[i];
        Vec3d p = body.pointAtEpoch(epoch);
        if(isnan(p.x)) continue;
        Draw3D::drawPointCross( (p-ref_pos)*zoom, sz );
        //glVertex3f((p.x-ref_pos.x)*zoom,(p.y-ref_pos.x)*zoom,(p.z-ref_pos.x)*zoom );
    }
}

void orbit( int n, const Orbit& orbit, double ang0, double ang1 ){
    Vec3d ps[n+1];
    double dang = (ang1-ang0)/n;
    orbit.toPoints( ang0, dang, n+1, ps );
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<=n; i++){
        glVertex3f((ps[i].x-ref_pos.x)*zoom,(ps[i].y-ref_pos.x)*zoom,(ps[i].z-ref_pos.x)*zoom );
    }
    glEnd();
}

void orbit_epochs( int n, const Orbit& orbit, double epoch0, double epoch1 ){
    Vec3d ps[n+1];
    orbit.toPoints_epochs( epoch0, epoch1, n+1, ps );
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<=n; i++){
        glVertex3f((ps[i].x-ref_pos.x)*zoom,(ps[i].y-ref_pos.x)*zoom,(ps[i].z-ref_pos.x)*zoom );
    }
    glEnd();
}

void boulder(const Boulder& b){
    //printf( "boulder span (%g,%g,%g)\n", b.span.x, b.span.y, b.span.z );
    glPushMatrix();
    Draw3D::rigidTransform( b.pos, b.rotMat, b.span );
    SphereSampling::drawIcosaMap( (Vec2i){b.nsamp,b.nsamp}, b.heights, b.hscale );
    glPopMatrix();
}

void rublePile(const RublePile& rb){
    for( const Boulder& b : rb.boulders){ boulder(b); }
}

// ================ From spaceTactics

inline Vec3d transformTrjPos( Vec3d* ps, int i){
    //printf("%f %f %f \n", ps[i].x,ps[i].y,ps[i].z );
    if( referenceBody ){
        if( bRefTrj ){ return (ps[i]-view_shift-referenceBody->trjPos[i])*view_scale; }
        else{          return (ps[i]-view_shift-referenceBody->pos      )*view_scale; }
    }else{             return (ps[i]-view_shift                         )*view_scale; }
}



void planet( const SpaceBody& b, int iTrj, double du ){
    Vec3d p;
    if(iTrj>0){
        p = b.getTrjPos(iTrj,du);
        if( referenceBody ) p.sub( referenceBody->getTrjPos(iTrj,du) );
    }else{
        p = b.pos;
        if( referenceBody ) p.sub(referenceBody->pos);
    }
    p.sub(view_shift);
    p.mul(view_scale);
    Draw3D::drawPointCross( p, 0.1 );
    Draw3D::drawSphere_oct(16, b.radius*view_scale, p );
}


void ship( const SpaceBody& b, int iTrj, double du ){
    Vec3d p;
    if(iTrj>0){
        p = b.getTrjPos(iTrj,du);
        if( referenceBody ) p.sub( referenceBody->getTrjPos(iTrj,du) );
    }else{
        p = b.pos;
        if( referenceBody ) p.sub(referenceBody->pos);
    }
    p.sub(view_shift);
    p.mul(view_scale);
    Mat3d mat = b.rotMat; mat.mul(b.sizes); mat.mul(scSz);
    Draw3D::drawMatInPos ( mat, p );
}


void trj( int n, const Vec3d * ps ){
    glBegin(GL_LINE_STRIP);
    Vec3f p;
    if( referenceBody ){
        if( bRefTrj ){ for(int i=0; i<n; i++){
            convert( (ps[i] - referenceBody->trjPos[i] - view_shift)*view_scale, p );
            //printf( "%i : %f %f %f \n", i, p.x, p.y, p.z );
            glVertex3d( p.x, p.y, p.z );
        }}else{ for(int i=0; i<n; i++){
            convert( (ps[i] - referenceBody->pos - view_shift)*view_scale, p );
            glVertex3d( p.x, p.y, p.z );
        }}
    }else{ for(int i=0; i<n; i++){
            convert( (ps[i] - view_shift)*view_scale, p );
            glVertex3d( p.x, p.y, p.z );
    }}
    glEnd();
}

void bodyTrj( const SpaceBody& b ){
    int n = trj_n;
    glBegin(GL_LINE_STRIP);
    Vec3d p,f;
    for(int i=0; i<n; i++){ p=transformTrjPos( b.trjPos, i ); glVertex3d( p.x, p.y, p.z ); };
    glEnd();
    if( b.trjThrust ){
        double sc = scF * view_scale * 1e+9;
        glBegin(GL_LINES );
        for(int i=iTrjMin; i<iTrjMax; i+=trjStep ){ p=transformTrjPos( b.trjPos, i ); f=b.trjThrust[i]; glVertex3d( p.x, p.y, p.z ); glVertex3d( p.x+f.x*sc, p.y+f.y*sc, p.z+f.z*sc );
            //if( &b == &world.ships[0] ) printf( "thrust %s[%i] (%f,%f,%f) \n", b.name.c_str(), i, f.x, f.y, f.z );
        };
        glEnd();
    }
}

void interactionTrj( const BodyInteraction& bi, const SpaceWorld& world ){
    int n = trj_n;
    glBegin(GL_LINE_STRIP);
    Vec3d p;
    glBegin(GL_LINES );
    for(int i=iTrjMin; i<iTrjMax; i+=trjStep ){
        p=transformTrjPos( world.ships[bi.i].trjPos, i ); glVertex3d( p.x, p.y, p.z );
        p=transformTrjPos( world.ships[bi.j].trjPos, i ); glVertex3d( p.x, p.y, p.z );
        //if( &b == &world.ships[0] ) printf( "thrust %s[%i] (%f,%f,%f) \n", b.name.c_str(), i, f.x, f.y, f.z );
    };
    glEnd();
}

} // namespace SpaceCrafting

#endif
