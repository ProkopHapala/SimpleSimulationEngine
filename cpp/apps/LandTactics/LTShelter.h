

#ifndef LTShelter_h
#define LTShelter_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"

#include "fastmath.h"
#include "tresholdFunctions.h"
#include "Vec2.h"
#include "Vec3.h"
#include "geom2D.h"
#include "Body2D.h"

#include "LTcommon.h"
//#include "LTUnitType.h"
//class LTFaction;

#include "LTUnit.h"


/*
TODO:
 - Shelters are points between static objects which are covered from particular directions
 - Unit can find sehlters automatically by sampling visibility lines in short-range neighborhood ( e.g. 10m..50m ? )
 - some units can enter perticular shelters ( e.g. infantery can enter house or bunker, tank can enter tank trench )


*/

enum LTSObjKind{
    tree, house, rocks
};

/*
class ShelterWall{ public:
    double armor;  // mm
    double filled; // percentage of filled area - how much holes there is
};
*/

class LTObjectType {  public:
    //Vec2d dims;            // lenght width
    LTSObjKind kind;
    int glo;
    //ShelterWall walls[4];
    //double height;        // above terrain?

    virtual void render( Vec2d pos, Vec2d dir ){};
    void render_glo( ){
        glo = glGenLists(1);
        glNewList(glo, GL_COMPILE);
            render( {0.0,0.0}, {1.0,0.0} );
        glEndList();
    };
};

class LTRectHouseType : public LTObjectType{
    Vec2d span = (Vec2d){10.0,5.0};

    virtual void render( Vec2d pos, Vec2d dir ){
        Vec2d dirT; dirT.set_perp(dir);
        dir.mul (span.a);
        dirT.mul(span.b);
        //printf( " (%f,%f) (%f,%f) (%f,%f) \n", pos.x, pos.y, dir.x, dir.y, dirT.x, dirT.y );
        glBegin(GL_LINE_LOOP);
            glVertex3f( pos.x-dir.x-dirT.x, pos.y-dir.y-dirT.y, 0.0 );
            glVertex3f( pos.x-dir.x+dirT.x, pos.y-dir.y+dirT.y, 0.0 );
            glVertex3f( pos.x+dir.x+dirT.x, pos.y+dir.y+dirT.y, 0.0 );
            glVertex3f( pos.x+dir.x-dirT.x, pos.y+dir.y-dirT.y, 0.0 );
        glEnd();
    };
};

class LTStaticObject{ public:
    // this object should represent Trees, Houses, Rocks ...
    //LTSObjKind kind;
    //int id; // void * ptr;
    Vec2d  pos    = (Vec2d){0.0,0.0};
    double radius = 10.0;
    Vec2d  dir    = (Vec2d){1.0,0.0};
    LTObjectType* type = NULL;
    //Vec3d span;

    void view(){
        Draw2D::drawCircle_d  ( pos, radius, 16, false );
        Draw2D::drawVecInPos_d( dir*radius, pos );
        Draw2D::drawShape( pos, dir, type->glo );
    };
};

/*
class LTShelter{ public:
    KinematicBody2D pose;  // 5 floats
    //Vec3d pose;       // pos angle (3 floats)
    LTShelterType* type     ;
    //uint8_t wallHP[4];   // wall health status
};
*/

#endif
