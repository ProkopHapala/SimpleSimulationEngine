

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

class ShelterWall{ public:
    double armor;  // mm
    double filled; // percentage of filled area - how much holes there is
};

class LTShelterType {  public:
    //Vec2d dims;            // lenght width
    ShelterWall walls[4];
    double height;        // above terrain?
};

class LTShelter{ public:
    KinematicBody2D pose;  // 5 floats
    //Vec3d pose;       // pos angle (3 floats)
    LTShelterType* type     ;
    //uint8_t wallHP[4];   // wall health status
};

#endif
