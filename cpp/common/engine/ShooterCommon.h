

#ifndef ShooterCommon_h
#define ShooterCommon_h

class ProjectileType{ public:
    double mass = 1.0;
    double explosiveMass = 0.0;
    double caliber       = 0.01; // [m]
    double crossection;
    double balisticCoef;

    inline void updateAux( double dragCoef ){
        crossection = M_PI*0.25*caliber*caliber;
        balisticCoef = 0.5*dragCoef*crossection/mass;
    }
    //double penetration0  = ;
};

#endif
