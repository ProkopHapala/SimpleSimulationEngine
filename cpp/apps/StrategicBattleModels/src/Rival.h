#ifndef Rival_h
#define Rival_h

enum AIR_ORDERS   {
    BELOW_CLOUDS = (1u << 0),
    MAXIMIZING = (1u << 1),
    MINIMIZED  = (1u << 2),
    MAXIMIZED  = (1u << 3),
};

enum GROUD_ORDERS {
    MINIMIZING = (1u << 0),
    MAXIMIZING = (1u << 1),
    MINIMIZED  = (1u << 2),
    MAXIMIZED  = (1u << 3),
};

struct UnitBrand{
    UnitType* type;
    int n;
};

struct GunType{
    float penetration; // [mm]
    float damage;      // [J?]
}

struct HitBox{
    float size;   // m
    float armor;  // armor
    float HPs;    // ammout of absorbed damage
    //float KillProb;
    
    float repairTime;
    float repairCost;
    
    //float vital;  // how much vital it is for the whole vehicle
}


inline float vehicleFunction(int n, float* componentHPs){
    // function of the whole vehicle requires function of all vital components
    // therefore it is products of function of sub components
    float function = 1.0;
    for(int i=0; i<n; i++}{ function *= componentHPs[i]; }
}


//inline float pentrationFunction( float pen){
//    return smoothstep();
//}


inline float hitToBoxes( int n, GunType& gun, HitBox* boxes, repairCost, repairTime ){
    float repairCost=0;
    float repairTime=0;
    float nkill     =0;
    for(int i=0; i<n; i++}{
        float dpen = boxes[i].armor - gun.penetration;
        //float npen = pentrationFunction(dpen/boxes[i].armor);
        float npen       = smoothstep( dpen/boxes[i].armor );
        float dmg        = npen/boxes[i].HPs;
        float killProb   = smoothstep( dmg );
        float nrepair    = dmg*(1-killProb);
        nkill      += killProb;
        repairCost += boxes[i].repairCost*nrepair;
        repairTime += boxes[i].repairTime*nrepair;
    }
    return nkill;
}



struct UnitType{
    float speed;
    float hitBoxes[];
    GunType guns[2];
}

struct AircraftUnitType{
    float speed;
    float climb;
    float ;
    GunType guns[2];
}




struct Rival{
    uint64_t order = 0;
    //std::vectro<UnitBrand> artilery; // this is rather battlerfield condition
    std::vectro<UnitBrand> support;
    std::vectro<UnitBrand> spearhead;
    std::vectro<UnitBrand> bombers;
    std::vectro<UnitBrand> fighters;
};

#endif
