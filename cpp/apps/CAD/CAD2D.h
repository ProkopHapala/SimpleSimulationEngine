
#ifndef  CAD2D_h
#define  CAD2D_h

//#include <vector>

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"
//#include "Mat3.h"
//#include "quaternion.h"
#include "geom2D.h"

namespace CAD{



class Blueprint{ public:
    std::vector<Ray2d   >  lines;   //
    std::vector<Circle2d>  circles; // circular args
    std::vector<Arc2d>     arcs;

    void addArc(const Vec2d& p1, const Vec2d& p2,  const Vec2d& p3 ){
        Circle2d c; c.from3points(p1,p2,p3);
        circles.push_back( c );
        Arc2d  arc;
        arc.fromCenter2points( c.p0, p1, p2 );
        arcs.push_back(arc);
    }

    void draw(){

    }

};



} // namespace SpaceCrafting

#endif
