

#ifndef LTUnit_h
#define LTUnit_h

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
#include "LTUnitType.h"
//#include "LTFaction.h"

const int Unit_JOB_mask_move   = 0b00010000;
const int Unit_JOB_mask_pos    = 0b00100000;
const int Unit_JOB_mask_LTUnit = 0b01000000;
const int Unit_JOB_mask_auto   = 0b10000000;

const int Unit_JOB_IDLE         = 0; // do whatever you want, regenerate at maximum speed
const int Unit_JOB_COVER        = 1; // duck, cover, hide, do not move

const int Unit_JOB_GOTO         = 2; // go to target with standard readyness, oppurtunistic shot at any target
const int Unit_JOB_ATTACK       = 3; // go to target while shooting; prefer cover, suppresive fire and oppurtinity to destroy enemy at target aread
const int Unit_JOB_HUNT         = 4; // chase LTUnit in order to kill it (move and fire)
const int Unit_JOB_SNEAK        = 5; // go to target carefully, and avoid any shooting
const int Unit_JOB_ASSAULT      = 6; // go to target at maximum speed; prefer movement over shooting

const int Unit_JOB_FIRE_AT_GOAL = 7;  // fire at position in terrain (do not move )
const int Unit_JOB_FIRE_AT_UNIT = 8;  // fire at particular LTUnit ( do not move )
const int Unit_JOB_FIRE_AT_ARK  = 9;  // fire at any enemy which appears in the target; choose priority tragets automatically
const int Unit_JOB_FIRE_AT_WILL = 10; // fire at any target in range;  choose priority tragets automatically

class LTFaction;

class LTUnit : public RigidBody2D {
	public:

    int n=100, ntot=100;
    double suppressed = 0.0; // between one and zero

    LTFaction  * faction;
    LTUnitType * type = NULL;

    Vec2d   goal;
    LTUnit  * opponent;
    //double    opponent_score = 0.0;

    Vec2d  attentionDir;
    double maxwf = 1.0;

    double radius = 15.0;

    int  job         = Unit_JOB_IDLE;
    int  default_job = Unit_JOB_IDLE;

    double reload = 0.0;
    int    ammo   = 100;

    // ===== inline functions

    void move_to_goal ( double dt );
    void fire_at_unit ( LTUnit * target );
    void update       ( double dt );
    double damage_ramp( double att, double def );
    void getShot    ( const Vec2d& from, int nburst, double area, double damage_const, double damage_kinetic, double penetration );
    void setOpponent( LTUnit * opponent_ );
    void setGoal    ( const Vec2d& goal_ );
    void setType    ( LTUnitType* type_ );
    void render     ( uint32_t c );
    void renderJob  ( uint32_t c );

    LTUnit(){};
    LTUnit( LTUnitType* type_, LTFaction* faction_, const Vec2d& pos_ );

};

#endif
