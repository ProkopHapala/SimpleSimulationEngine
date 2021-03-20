
#ifndef LTrender_h
#define LTrender_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"

#include "LTcommon.h"
#include "LTUnitType.h"
#include "LTUnit.h"

void render( const LTUnit& u, uint32_t color, int iLOD ){
    Draw::setRGBA(color);
    //printf( "unit.type->kind %i %s  (%f,%f)  (%f,%f) \n", type->kind, sUnitKind[type->kind], pos.x,pos.y,  rot.x,rot.y );
    switch(u.type->kind){
        case(UnitKind::inf):
            Draw2D::drawRotT      (u.pos, u.rot, {u.type->sz.a, u.type->sz.b} );
            break;
        case(UnitKind::gun):
            Draw2D::drawRotTriangle(u.pos, u.rot, {u.type->sz.a*0.5, u.type->sz.b});
            Draw2D::drawVecInPos_d (u.turret_dir*u.type->sz.a,u.pos);
            break;
        case(UnitKind::tank):
            Draw2D::drawRotRect   (u.pos, u.rot, {u.type->sz.a, u.type->sz.b});
            Draw2D::drawCircle_d  (u.pos, u.type->sz.b*0.5,16,false);
            Draw2D::drawVecInPos_d(u.turret_dir*u.type->sz.a,u.pos);
            break;
        case(UnitKind::stug):
            Draw2D::drawRotRect   (u.pos, u.rot, {u.type->sz.a, u.type->sz.b});
            Draw2D::drawVecInPos_d(u.rot*u.type->sz.a,u.pos);
            break;
    }
}


void render( const LTSquad& sq, uint32_t color, int iLOD, bool bDrawGoal ){
    //printf( "squad \n" );
    //printf( "squad pos (%f,%f) \n", pos.x, pos.y );
    //glColor3f( c.x, c.y, c.z );
    Draw::setRGBA(color);
    Draw2D::drawCircle_d( sq.pos, sq.radius, 16, false );
    Draw2D::drawVecInPos_d( sq.attentionDir*sq.radius*0.5, sq.pos );
    Draw2D::drawVecInPos_d( sq.rot*sq.radius,              sq.pos );
    //printf( " %f %f \n", attentionDir.x, attentionDir.y );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
    if(iLOD>0){
        for(const LTUnit& u: sq.units ){
            //printf( "unit \n"  );
            render( u, color|0xFF000000, iLOD );
            Draw::setRGBA(color);
            if( bDrawGoal ){
               Draw2D::drawLine_d( u.goal_pos, u.pos );
            }
        }
    }

    char str[8];
    sprintf(str,"%4i",sq.n);
    //Draw2D::drawString( str, (float)pos.x, (float)pos.y, 0.4f, default_font_texture );
    Draw2D::drawText( str, 0, {sq.pos.x,sq.pos.y}, 0.0, default_font_texture, 2.0 );
}

void renderJob( const LTSquad& sq, uint32_t c){
    if(sq.job == Unit_JOB_GOTO         ) Draw2D::drawLine_d( sq.pos, sq.goal );
    if(sq.job == Unit_JOB_FIRE_AT_UNIT ) Draw2D::drawLine_d( sq.pos, sq.opponent->pos );
    //printf( "render (%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", pos.x, pos.y, pos.x, c.x, c.y, c.z );
}

#endif
