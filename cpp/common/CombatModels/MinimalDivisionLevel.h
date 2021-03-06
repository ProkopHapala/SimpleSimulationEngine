
/*

ToDo:
 - Moral system (organization)
 - Advance system
 - Limit concentration
    - Higher losses from artilery when concentrated
 - Surprise - firepower of defender in first round is limited
 - Limit Range by Terrain ( in city or forrest cannot shoot long distances )
 - Fire-Rate vs FirePower
   - Some weapons can shoot lot of projecties with low precission(SMG), while other shoot few projectiles precisely (sniper)
       - should consider in ranged combat

=======================

Model is based on ratio of several different kinds of units:
- Infantry
- Tanks (Armor)
- Support:
 - Artilery
 - Aircraft

## Units ability can be modified by equipment
 - Transports  - increase transport speed of Infantery and Guns
 - AT Guns     - imporve anti-armor capability of infantery, decrease speed, increase transport weight and supply consumption
 - Machineguns - increase light firepower, 
 - Mortars    -

Armor:    none(Infantery, Tucks), light(APC,), heavy( tanks, bunkers )
Mobility: none(bunkers), slow( Infantery, Artilery ), fast( Artilery, SPG )

## Parameters Of Units:
 - Supply consumption - strain of supply chain and production [tons/day] (may distiguis Fuel, Ammo, Food )
 - Transport Weight   - how much transport capacity [ton-day/km] used when transfered by infrastructure
 - Combat Movement Speed  - Movement Speed when deployed for fight 
 - Transit Movement Speed - Movement Speed when undeployed
 - Firepower
   - Light(Soft) - against un-armored targets
   - Medium      - Against APC and AFV
   - Heavy       - Against Tanks
   - SuperHeavy  - Bunkers and Battle ships 
 - CQC ( Close-Quarters Combat bonus onf Firebpower )
    - FirePower
       - Light
       - Medium
       - Heavy
       - Super Heavy
 - Cover bonus
   - 

## Parameters of Battlefield / Battle modifiers
 - cover - typically vavours light units like infantery
 - Recon - 

# Fortification
 - Directional - each side of cell (typicaly Hex Cell) can be fortified, The fortification is directional and not transitive - fortification of cell A vs B does not provide any bonus of B vs A

## Waves/Phases of battle (Prvo-sledove / Druho-sledove oddily)
 - First line tipically suffer more heavy losses, therefore units which are more valuable (either more costly, or required elsewhere) should come second
 - First line also cause enemy to commit and engage, while second line deal the final blow
 - fist line should be tipically 
 - Two schemes
    - Armor First
        - Sufficient Armor attack is able to penetrate enemy line very quickly, without actually eliminating all defenders
        - Armor provide considerable suppresive fire which cover advance of infantery
        - One Enemy Line is penetrate infantery can advance much more easily and clean up the remaining defenders, and capturing lot of prissoners of war (POWs)
        - However, if the armor is insufficient to penetrate, e.g. if enemy posses enough anti-tank weapons, the whole armore wave can be eliminated which is extremely costly
        - If the enemy is well prepared, and our recon is insufficient, than blind armored atack most probalbly fails, because armored units does not have enough recon capability and enoguh time locate and eliminate AT weapons
    - Infantery First
        - Cautious attack of infantery 

# Damage Model Phases
 - Assemble Firepower Potential - weighted sum of firepower in each cathegory (light,heavy ... )
 - Assemble Visibility - total visibility of units in each catheogory, calculate share assigned to each units (normalized by total visibility)
    - if visibility is under certain level the enemy reaction is not sparked. Firepower potential is not fully used in that case, since no target is avilable.
 - Evaluate Effective Firepower (i.e. that which is really used in combat)
 - Distribute Firepower among target groups (different kinds of units on the battle field) according to visibility share
 - Actuall damage model - firepower has diffenrent effect on targets depending on
   - hit probability         - depending on size and cover (strongly related but not nessesarily the same as visibility; e.g. fortification decrease hit probability but not visibility)
   - penetration probability - depending on armor and armor penetration rating of weapons (AP)
   - damage susceptibiility / after penetration damage (APD)  - depend on both weapon and target. There is often compromise between armor penetration (AP) and after penetration damage (APD)
   - Effect : levels Hard-kill, soft-kill
        - Suppressed  - temporarily does not contribute to the battle
        - Demoralized - decrease combat activity, increas probability to surrender
        - Disabled    - permanently eliminated from the battle until repaired/healed
        - Killed      - completely eliminated from list of forces

# Damage distribution
 - Fire-power is distributed according to visibility - i.e. most visible units (lime tanks) attract most fire
 - Units which can provide both light and heavy fire distribute have to choose (it is exclusive), which is done also depending on visibility

## Concept of granularity
* Cost of units must be considered

## Ammo Saving mode
 * Consumption of ammo can be HUGE, it make sense to balance how much confidence to hit target justify to fire the weapon


* UnitTypes
Inf
Inf_MG
Inf_AT
Inf_ATGM (MANPATS)
Inf_AAGM (MANPADS)

#### Infantery Modification

MG:
 - huge firepower against infantry in defence
 - significant supperssion
 - not as efficient in offence
 - slightly decreased movement speed
 - huge ammo consumption

Mortar:
 - serious increase of indirect suppression power (cheap artilery)
 - increase firepower against fortification
 - way less efficient in direct fire than machine gun or even rifle, slow movement
 - huge ammo consumption

Sappers/FlameTrower:
 - huge increase of firepower in close-quarter combat
 - only CQC; Limited range, useless at long range

AT:
 - reasonable firepower against lighter tanks in defence
 - increase firepower against fortified infantery
 - seriously hampered combat_speed

RPG:
 - considerably improved Anti-Tank and anti-fortification capability (splash damage)
 - only CQC

ATM:
 - considerable Anti tank capability
 - costly - single-shot (no ammo re-charge)

AAM:
 - considerable anti-helicopter capability

Transport ()
 - increase transit speed

Armored Transport (APC)
 - increase anti-personal capability (heavy MG with large ammo supply)
 - increase combat speed
 - worse cover/stealth

Armored Transport (AFV)
 - increase anti-tank and anti-personal capability
 - increase combat speed
 - worse cover/stealth
 - costly

#### Armored Units

Tank
 - heavy armor (can be killed only by very limited range of weapons)
 - heavy firepower, especially against other heavy
 - good combat speed
 - costly

Hell
 - huge anti-tank capability
 - huge combat speed
 - ignores terrain and fortification
 - ignores artilery
 - easily killed by AAM - no-repair
 - costly

#### Support units

Artilery
 - Provide indirect fire from long range (30km)
 - significantly help suppress enemy, break mass enemy attack, and siege fortification
 - Two modes:
    - Point target - target have to be marked by ground troops
    - Area - low hit probability unless enemy is concentrated
 - no direct combat capability.
 - inefficient vs armor
 - inefficient vs. covered units
 - extreme ammo consumption when used recklessly

Close Air Support
 - Can be quickly concentrated from far distance
 - Modern versions Can pracically any target
 - Target had to be marked by ground units first
 - huge fuel consumption when used recklessly
 - very costly - especially when destroyed

*/

#ifndef MinimalDivisionLevel_h
#define MinimalDivisionLevel_h

int iDEBUG=0;


inline double firePowerSaturation( double fp, double vis ){
    //return vis*erf(fp/vis);
    if(fp==0) return fp;
    return erf(vis/fp);
    //  firepower cannot be higher than visibility
}

inline double supress_func(double shots){
    return exp(-shots);
}

inline double logistic(double x){
    return 1/( 1 + exp(-x) );
}



/*
inline double fpDecay( double  dist, double reach ){
    double  c = dist/reach;
    return  1 - c*c;
}
*/

// ==================================
// =====   BattleField
// ==================================

struct BattleField{
    // weather
    //  - air - cloudy
    //  - fog
    //  - mud
    // terrain
    //double viewRate  = 5;
    double viewRange = 10000; // [m]
    double maxCamo   = 0;     // [0..1]
};

// ==================================
// =====  WeaponType
// ==================================

struct WeaponType{
    double direct_reach;    // distance
    double indirect_reach;       // 
    double indirect_supression;  // supression in idirect or direct fire
    double movePenalty;  // [0.0..1.] decreased firepower when moving
    double supply;       // supply consumed when fireing [/attack round]
    double firepower[4]; // against armor-class : light[<10mm], medium[<50mm], heavy[<500mm], super heavy[>500mm]
    //double damage   [4]; // damage caused to diven armor-class

    int fromString(const char* s){
        return sscanf( s, "%lf %lf %lf %lf %lf    %lf %lf %lf %lf  \n",
            &direct_reach, &indirect_reach, &indirect_supression, &movePenalty, &supply,
            &firepower[0],&firepower[1],&firepower[2],&firepower[3]
            //&damage[0],   &damage[1],   &damage[2],   &damage[3] 
        );
    }

    int info(char* s){
        return sprintf( s, "reach %g / %g %g \n movePen %g \n supply %g \n fire[%g,%g,%g,%g] \n ", 
            direct_reach, indirect_reach, indirect_supression, movePenalty, supply,
            firepower[0],firepower[1],firepower[2],firepower[3]
            //damage[0],   damage[1],   damage[2],   damage[3] 
        );
    }

    inline double getDirectDecay(double dist){
        double  c = dist/direct_reach;
        if(c>1) return 0;
        return  1 - c*c;
    }
};

// ==================================
// =====  UnitType
// ==================================

struct UnitType{
    int    armorClass;     // light[<10mm], medium[<50mm], heavy[<500mm], super heavy[>500mm]
    double size;           // [m]   => visibility and hit probability
    double supply;         // [kg/day]
    double weight;         // strategic   speed
    double baseCost;
    double speed_transit;  // [km/h] operational speed
    double speed_combat;   // [km/h] tactical    speed
    WeaponType* primary;    // this weapon is prefered
    WeaponType* secondary;  // used only when primary cannot be used

    int fromString(const char* s){
        return  sscanf( s, "%i %lf %lf %lf %lf   %lf %lf  \n", 
            &armorClass,
            &size, &weight, &supply, &baseCost, 
            &speed_transit, &speed_combat   
            //primary, secondary,
        );
    }

    int info(char* s, bool bDeep){
        char *s0 = s;
        s+=sprintf( s, " armor %i size %g weight %g \n supply %g cost %g \n speed %g / %g \n", 
            armorClass,
            size, weight, supply, baseCost, 
            speed_transit, speed_combat   
            //primary, secondary,
        );
        //printf( " s-s0 %i \n", s-s0);
        if(bDeep){
            if(primary  ){ s+=sprintf(s,"primary:\n"  );   s+=primary  ->info(s); }
            if(secondary){ s+=sprintf(s,"secondary:\n");   s+=secondary->info(s); }
        }
        return s-s0;
    }

    inline double battlefield_speed(const BattleField& cond){
        return speed_combat; // ToDo : effect of terrain and weather
    }

};

// ==================================
// =====  UnitState
// ==================================

struct UnitState{
    UnitType* type;
    int n;    // number combat-capable units
    int ntot; // total number of units including disabled
    double cost;
    // state / temp / aux
    double shot_attaction; // attractivity for enemy to target this unit under current tactical conditions
    double supressed;      // [0.0..1.0] is the unit suppressed by enemy fire ?
    double camo;           // [0.0..1.0] camouflage level ... how much more difficult it is to be seen
    double dug;            // [0.0..1.0] entranchement level ... how much more difficult it is to be hit (basically multiply size)
    double moral;          // [0.0..1.0]
    // orders
    double stealth;        // how much the unit risk to be uncovered, high-profile action e.g. run or fire decrease stealth
    double safeAmmo;       // [0.0..1.0] determine trade-off between ammo supply consumtion and firepower
    //double debug_got_shot;
    double got_fire;
    UnitState()=default;

void init( UnitType* type_, int n_ ){
    type=type_; n=n_; ntot=n;
    cost = type->baseCost;      // how valuable the unit is in curent tactical situation
    // state
    supressed  = 1.0; // [0..1] when under enemy fire, the movement and firepower is limited by this factor ( supressed=1.0 - no effect, supressed=0.0 full effect )
    dug        = 1.0; // [0..1] unit is hidden
    camo       = 1.0; //
    // orders
    safeAmmo   = 1.0;
    stealth    = 1.0; //how much stealthy the unit tries to be
}

UnitState(UnitType* type_, int n_){
    init(type_, n_); 
}

inline void reset(){
    moral     = 1;
    supressed = 1;
    n         = ntot;
}

inline double getFireActivity( ){
    return fmin( fmin( safeAmmo, supressed), moral*(2-moral) );
}

inline void takeFire( double firepower_got ){
    double invn = 1./n;
    got_fire    = firepower_got*invn; // ToDo - perhaps += in odrer to consider indirect fire and bombardment
    supressed   = supress_func(got_fire);
    //if(iDEBUG>0) printf( "firepower_got %g got_fire %g supressed %g \n", firepower_got, got_fire, supressed  );
}

inline void  updateDamage( const BattleField& cond, double dt){
    double  shots  = got_fire*dt;
    double dmoral = shots*0.05;
    int    dn     = n*shots*0.002;
    moral   -= dmoral; shots*0.05;   if(moral<0)moral=0;  // ToDo: replace magic-factor  0.05 by something more justified
    n       -= dn;                   if(    n<0)n=0;
    if(iDEBUG>0) printf( "updateDamage dt %g got_fire %g dmoral %g dn %i  \n", dt, got_fire, dmoral, dn );
}

inline void recover( const BattleField& cond, double dt ){
    moral += (1.1-moral)*0.08*dt; if(moral>1)moral=1;
}

inline double advance_speed( const BattleField& cond ){
    double suppress2 = supressed*supressed;
    //printf( "suppress2 %g \n", suppress2 );
    return type->battlefield_speed(cond)*suppress2;
}

inline double visibilityFunction( double dist, double zoom, const BattleField& cond ){
    // maxCamo determined by terrain
    double invAng   = dist/(zoom*type->size);
    //double decay    = 1/(1+invAng*invAng);
    //double fcut   = 1/(1+exp( cond.viewRate*(dist-cond.viewRange)/cond.viewRange) );  // ToDo : This can be factore out
    double argCut   = dist/cond.viewRange;
    //double fcut     = 1/(1+argCut*argCut ); 
    double decay    = 1/((1+invAng*invAng)*(1+argCut*argCut ));
    //double hidden = 1. - maxCamo*camo*stealth;
    double hidden = cond.maxCamo*camo*stealth;
    //printf("visibilityFunction: cross %g(size %g invDist2 %g) * hidden %g (maxCamo %g camo %g stealth %g) \n",  cross, type->size, invDist2, hidden, maxCamo, camo, stealth );
    return decay*hidden;
    //return fcut;
    //return decay*hidden;
}

};

// ==================================
// =====   Division
// ==================================

struct Division{
    int n_inf; // number of infantery
    UnitType transport;
    std::vector<UnitState*> units;

    void reset(){
        for(UnitState* unit : units){ unit->reset(); }
    }

};



// ==================================
// =====   CombatSide
// ==================================

struct CombatSideTemp{
    double ammo_spend;
    // prec shots / vs suppess/burst/scattered shots
    double primary_firepower[4];
    double secondary_firepower[4];
};

struct CombatSide{
    Division composition;
    double zoom            = 50; // zoom enemy unit size in visibility calculation (i.e. better target search on distance)
    double firepower        [4]; // total amount of firepower in each category {light, medium, heavy, superheavy}
    double target_attraction[4]; // total visibility in each category          {light, medium, heavy, superheavy}
    double totalVisibility;
    double advanced;             // how much the unit advanced or retreated form target line
    double target_dist = 0;      // distance to which we try to advance

void accumUnit(UnitState& unit, double dist, double zoom, const BattleField& conds, CombatSideTemp& tmp ){
    // ---- shot Attraction
    //double invDist2 = 1/(dist*dist);
    double visibility = unit.visibilityFunction( dist, zoom, conds );
    //printf( "visibility %g \n", visibility );
    totalVisibility  += visibility;
    //double attaction  = attractivityFunction( unit.type.size, unit.value, unit.camo, unit.stealth );
    double hitProb      = erf( visibility );
    //printf( "visibility %g hitProb  %g \n", visibility, hitProb  );
    unit.shot_attaction = unit.cost * hitProb;
    double shot_needed  = unit.n * unit.shot_attaction;
    target_attraction[ unit.type->armorClass ] += shot_needed;
    //printf( "visibiliy %g hitProb  %g shot_needed %g n %i cost %g \n", visibility, hitProb, shot_needed, unit.n, unit.cost );
    //if(iDEBUG>0) printf( "armor[%i] visibility %g hitProb %g cost %g attraction %g \n" , unit.type->armorClass, visibility, hitProb, unit.cost, unit.shot_attaction );
    //return;

    // ---- FirePower
    //const UnitState& unit  = *unit_;
    double activity  = unit.n * unit.getFireActivity();
    double safeAmmo2 = unit.safeAmmo*unit.safeAmmo;
    
    const double *fps1,*fps2;
    double fdec1,fdec2; 
    fps1  = unit.type->primary->firepower;
    fdec1 = unit.type->primary->getDirectDecay( dist );
    //const double* fps2;
    if(unit.type->secondary){ 
        //fps2 = unit.type->secondary->firepower * fpDecay( dist, unit.type->secondary->direct_reach );
        fps2  = unit.type->secondary->firepower;
        fdec2 = unit.type->secondary->getDirectDecay( dist );;
    }else{ fps2 = 0 ;};
    //printf( " activity %g fp1,2 %g,%g  \n", activity, fps1, fps2 );
    for(int i=0; i<4; i++){
        double fp1,fp2;
        double dfp1,dfp2;
        fp1  = fps1[i];
        dfp1 = activity * fp1 * fdec1;
        tmp.primary_firepower[i] += dfp1;
        if(fps2){
            fp2  = fps2[i];
            dfp2 = activity * fp2 * fdec2;
            tmp.secondary_firepower[i] += dfp2;
        }else{
            fp2=0;
        }
        tmp.ammo_spend += (fp1 + fp2)*safeAmmo2;
        //if(iDEBUG>0) printf( "->fp[%i] activ %g fp(%g,%g) fdec(%g,%g) dfp(%g,%g) \n", i, activity,  fp1,fp2,  fdec1,fdec2,  dfp1,dfp2 );
    }
}

// Maybe call it rather Fire-Attractivity or Fire-Magnet
void gatherFire( const BattleField& conds, double dist, double zoom, CombatSideTemp& tmp ){
    //double invDist2 = 1/(dist*dist);
    //totalVisibility=0;
    tmp.ammo_spend =0;
    for(int i=0; i<4; i++){ 
        tmp.primary_firepower  [i]=0;
        tmp.secondary_firepower[i]=0;
        target_attraction      [i]=0; 
    }
    for(UnitState* unit : composition.units){
        //printf("unit->type->armorClass %i \n", unit->type->armorClass );
        accumUnit( *unit, dist, zoom, conds, tmp );
    }
}

inline void match( const double* attraction, CombatSideTemp& tmp ){
    for(int i=0; i<4; i++){
        double fp1 = tmp.primary_firepower  [i];
        double fp2 = tmp.secondary_firepower[i];
        double vis = attraction[i];
        double f1  = firePowerSaturation( fp1,     vis );
        double f2  = firePowerSaturation( fp1+fp2, vis );
        firepower[i] = f1*fp1 + f2*fp2;
        //if(iDEBUG>0) printf( "->match[%i] vis %g fp(%g,%g) f(%g,%g) dfp(%g,%g) \n", i, vis, fp1, fp2, f1, f2, fp1*f1, fp2*f2 );
    }
}

/*
void distributeFirePower( const double* firepower, double dt ){ 
    double normalized_firepower[4];
    for(int i=0; i<4; i++){
        double att = target_attraction[i];
        if(att>0){ normalized_firepower[i] = firepower[i] / att; }else{ normalized_firepower[i]=0; };
        //normalized_firepower[i] = firepower[i] * invAttraction[i]; 
        //if(iDEBUG>0) printf( "normalized_firepower[%i] %g <= fp %g att %g \n", i, normalized_firepower[i], firepower[i], att );
    };
    for(UnitState* unit : composition.units){
        int iclass = unit->type->armorClass;
        double firepower_got = unit->n * unit->shot_attaction * normalized_firepower[iclass];
        unit->getFire( firepower_got, dt ); // internal damage model
        //if(iDEBUG>0) printf( "unit.armor[%i] got_fp %g attraction %g norm_fp %g \n", iclass, firepower_got,  unit->shot_attaction, normalized_firepower[iclass] );
    }
}

double tryAdvance( const BattleField& conds, double dt ){  // ToDo: this may be possible to integrate within distributeFirePower 
    double vmin = 1e+300;
    for(UnitState* unit : composition.units ){
        double vi = unit->advance_speed(conds);
        vmin = fmin(vmin,vi);
    }
    printf( "speed %g dt %g \n", vmin, dt );
    // ToDo: if the advance is hampered by single exceptionally slow unit, it may be reasonable to assign it to second wave ????
    return vmin*dt;
}
*/

double moveUnderFire( const BattleField& conds, const double* firepower, double dt_max, double dist ){ 
    double normalized_firepower[4];
    for(int i=0; i<4; i++){
        double att = target_attraction[i];
        if(att>0){ normalized_firepower[i] = firepower[i] / att; }else{ normalized_firepower[i]=0; };
        //normalized_firepower[i] = firepower[i] * invAttraction[i]; 
        //if(iDEBUG>0) printf( "normalized_firepower[%i] %g <= fp %g att %g \n", i, normalized_firepower[i], firepower[i], att );
    };
    double vmin = 1e+300;
    //printf("----- %i \n", composition.units.size() );
    for(UnitState* unit : composition.units){
        int iclass = unit->type->armorClass;
        double firepower_got = unit->n * unit->shot_attaction * normalized_firepower[iclass];
        //if(iDEBUG>0) printf( "unit->n %i ", unit->n );
        unit->takeFire( firepower_got ); // internal damage model
        double vi = unit->advance_speed(conds);
        vmin = fmin(vmin,vi);
        //if(iDEBUG>0) printf( "unit.armor[%i] got_fp %g attraction %g norm_fp %g \n", iclass, firepower_got,  unit->shot_attaction, normalized_firepower[iclass] );
    }
    double dt_need = (dist-target_dist)/vmin; // time the wave need to advance to target distance
    double dt      = fmin( dt_need, dt_max );
    advanced+=vmin*dt;
    //if(iDEBUG>0) printf( "speed %g dt %g ddist %g advanced %g \n", vmin, dt, dt*vmin, advanced );
    return dt;
}

void updateDamage( const BattleField& conds, double dt ){ 
    for(UnitState* unit : composition.units){  // this can be esily moved under "Division" class
        unit->recover     ( conds, dt ); 
        unit->updateDamage( conds, dt );
    }
}

void clearTemp(){ 
    for(UnitState* unit : composition.units){
        unit->supressed = 1;
        unit->moral     = 1;
    }
}

};

// ==================================
// =====   CombatSide
// ==================================

struct Combat{
    double dist;
    BattleField conds;
    CombatSide attacker;
    CombatSide defender;

double round(double dt_max, double advance_dist ){
    /*
    attacker.assembleTargetAttraction( conds, dist );
    defender.assembleTargetAttraction( conds, dist );
    attacker.assembleFirePower( defender.target_attraction );
    defender.assembleFirePower( attacker.target_attraction );
    */
    //double attfp1[4],attfp2[4],deffp1[4],deffp1[4];
    //if(iDEBUG>0)printf( "=======combat round dist %g \n", dist );
    CombatSideTemp atttmp, deftmp;
    //double invDist2 = 1/(dist*dist);
    //if(iDEBUG>0)printf( "attacker::gather \n");
    attacker.gatherFire( conds, dist, defender.zoom, atttmp );
    //if(iDEBUG>0)printf( "defender::gather \n");
    defender.gatherFire( conds, dist, attacker.zoom, deftmp );
    //return;
    //if(iDEBUG>0)printf( "attacker.match -> defender \n");
    attacker.match( defender.target_attraction,  atttmp );
    //if(iDEBUG>0)printf( "defender.match -> attacker \n");
    defender.match( attacker.target_attraction,  deftmp );

     

    // ---- distribute firepower and take damage
    //if(iDEBUG>0)printf( "attacker take fire \n");
    if(iDEBUG>0) printf("att \n"); 
    double dt = attacker.moveUnderFire( conds, defender.firepower, dt_max, advance_dist );
    //if(iDEBUG>0)printf( "defender take fire \n");
    if(iDEBUG>0) printf("def \n"); 
    defender.moveUnderFire( conds, attacker.firepower, 0, 0 );

    attacker.updateDamage( conds, dt );
    defender.updateDamage( conds, dt );

    dist = -(defender.advanced + attacker.advanced);
    return dt;
    // ---- movement
    //double advance_dist = attacker.tryAdvance(conds, dt ); // if remain some action points
    //printf( "advance_dist %g \n", advance_dist );
    //defender.tryRetreat(conds); // if defender is superior in ranged weapons it may be advantagenous to retreat, also if defenders moral break, it is forced to retreat
    //double dist = fmin( attacker.target_dist, dist-advance_dist );
}

double run(double dt_max, double fdist ){
    const int nMaxIter = 10;  
    for(int i=0; i<nMaxIter; i++){
        //attacker. fdist;
        //if(iDEBUG>0) printf( "combat.run[%i] dt_max %g dist %g advance_dist %g \n", i, dt_max, dist, fdist*dist );
        dt_max -= round( dt_max, fdist*dist );
        //if(iDEBUG>0) printf( "combat.run[%i] t_left %g dist %g \n\n", i, dt_max, dist );
        if(dt_max < 1.) break;                // time out
        if( attacker.advanced > 0 ) break;    // attacker reached the target
        //if( attacker.moral > 0 ) break;     // attacket retreat
        //if( defender.advanced > 0 ) break;  // defender retreat / surrender
    }
    return dt_max;
};

void start(double dist_){
    dist=dist_;

    attacker.clearTemp();
    defender.clearTemp();

    attacker.advanced = -dist;
    defender.advanced = 0;
}

void clear(){
    attacker.composition.units.clear();
    defender.composition.units.clear();
}

/*
// Maybe call it rather Fire-Attractivity or Fire-Magnet
double assembleTargetAttraction( const BattleField& conds, double dist ){
    double totalVisibility = 0;
    double invDist2 = 1/(dist*dist);
    for(UnitState* unit : composition.units){
        double visibility = unit->visibilityFunction( invDist2, conds.maxCamo );
        totalVisibility += visibility;
        //double attaction  = attractivityFunction( unit.type.size, unit.value, unit.camo, unit.stealth );
        unit->shot_attaction = unit->cost * erf( visibility );
        target_attraction[ unit->type->armorClass ] += unit->n * unit->shot_attaction;
        printf( "armor[%i] visibility %g cost %g attraction %g \n" , unit->type->armorClass, visibility, unit->cost, unit->shot_attaction );
    }
    return totalVisibility;
    // evaluate how attractive is each unit for enemy to shoot at
}

double assembleFirePower( const double* attraction ){
    double primary_firepower  [4];
    double secondary_firepower[4];
    double ammo_spend = 0;
    for(const UnitState* unit : composition.units){
        //const UnitState& unit  = *unit_;
        double activity  = unit->n * unit->safeAmmo * unit->supressed;
        double safeAmmo2 = unit->safeAmmo*unit->safeAmmo;
        const double* fps1 = unit->type->primary  ->firepower;
        const double* fps2 = unit->type->secondary->firepower;
        for(int i=0; i<4; i++){
            double fp1 = fps1[i];
            double fp2 = fps2[i];
            ammo_spend += (fp1 + fp2)*safeAmmo2;
            primary_firepower  [i] += activity * fp1;
            secondary_firepower[i] += activity * fp2;
        }
    }
    for(int i=0; i<4; i++){
        double fp1 = primary_firepower[i];
        double fp2 = secondary_firepower[i];
        double vis = attraction[i];
        fp1 *= firePowerSaturation( fp1,     vis );
        fp2 *= firePowerSaturation( fp1+fp2, vis );
        firepower[i] = fp1 + fp2;
    }
    return ammo_spend;
    // choose targets accorgint to attraction
    // if there is not enought visible/attractive targets, we try use secondary weapon
}
*/


};

#endif