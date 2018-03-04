

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

#include "geom2D.h"


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

class LTLinearObject{ public:
    int kind;
    Vec2d p1,p2;
    double width = 5.0;

    // TODO:
    // more sophisticated system of cover and terrain properties
    double cover=1.0; // TODO: maybe in future this shoudl be in Type ?

    inline Vec2d dp( const Vec2d& p ){ // derivative of distance squared
        /*
        double r2,r2min;
        Vec2d  dp,dpmin;
        dpmin=p-p1; r2min=dpmin.norm2();
        // linar segment
        Vec2d  d = p2-p1;
        double c = d.dot(dpmin)/d.norm(); // chan this be optimized ?
        if(c>0){
            c/=d.norm();
            dpmin.add_mul(d,-c);
            r2min = dpmin.norm2();
        }
        dp=p-p2;  r2=dp.norm2();  if(r2<r2min){ r2min=r2; dpmin=dp; }
        return dpmin;
        */
        return dpLineSegment(p,p1,p2);
    }

    inline char intersection( const Vec2d& op1, const Vec2d& op2, Vec2d& p ){
        //intersection_st( p1, p2, op1, op2, s, t );
        return intersection_point(  p1, p2, op1, op2, p );
    }
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
