

#ifndef PolyLineFormation_h
#define PolyLineFormation_h

#include <vector>
#include "SoftPolyLine2D.h"

class PolyLineFormation : public SoftPolyLine2D { public:
    // soldiers are in boxes attached to each line-segment
    // - they are attracted to the line segment, but otherwise can move a bit freely
    // - each line segment has its onw bounding box, which is spherical and lager than line segment
    // - there is also global rectanguler bounding box for the whole formation
    // - if soldier does not fit to particular line segment (e.g. if it is full, or of soldier is too far), he can be assigned to the global bbox
    //   - neighboring segments automatically exchange soldiers to balance themself
};

#endif
