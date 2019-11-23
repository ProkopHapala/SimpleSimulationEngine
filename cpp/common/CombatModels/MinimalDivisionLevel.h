
/*

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


double firePowerSaturation( double fp, double vis ){
    return erf(vis/fp);
    //  firepower cannot be higher than visibility
}


// ==============================================

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
};

// ==============================================

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

};

// ==============================================

struct UnitState{
    UnitType* type;
    int n;    // number combat-capable units
    int ntot; // total number of units including disabled
    double cost;
    double shot_attaction; // attractivity for enemy to target this unit under current tactical conditions
    double supressed;      // [0.0..1.0] is the unit suppressed by enemy fire ?
    double camo;           // [0.0..1.0] camouflage level ... how much more difficult it is to be seen
    double dug;            // [0.0..1.0] entranchement level ... how much more difficult it is to be hit (basically multiply size)
    double stealth;        // how much the unit risk to be uncovered, high-profile action e.g. run or fire decrease stealth
    double safeAmmo;       // [0.0..1.0] determine trade-off between ammo supply consumtion and firepower


void getFire( double firepower_got ){} 

double visibilityFunction( double invDist2, double maxCamo ){
    // maxCamo determined by terrain
    return (1. - maxCamo*camo*stealth )*type->size*invDist2;
}

};

// ==============================================

struct Division{
    int n_inf; // number of infantery
    UnitType transport;
    std::vector<UnitState> units;
};

// ==============================================

struct BattleField{
    // weather
    //  - air - cloudy
    //  - fog
    //  - mud
    // terrain
    double maxCamo;
};


// ==============================================

struct CombatSide{
    double firepower        [4]; // total amount of firepower in each category {light, medium, heavy, superheavy}
    double target_attraction[4]; // total visibility in each category          {light, medium, heavy, superheavy}
    Division composition;

// Maybe call it rather Fire-Attractivity or Fire-Magnet
double assembleTargetAttraction( const BattleField& conds, double dist ){
    double totalVisibility = 0;
    double invDist2 = 1/(dist*dist);
    for(UnitState& unit : composition.units){
        double visibility = unit.visibilityFunction( invDist2, conds.maxCamo );
        //double attaction  = attractivityFunction( unit.type.size, unit.value, unit.camo, unit.stealth );
        unit.shot_attaction = unit.cost * erf( visibility );
        for(int i=0; i<4; i++) target_attraction[ unit.type->armorClass ] += unit.n * unit.shot_attaction;
    }
    return totalVisibility;
    // evaluate how attractive is each unit for enemy to shoot at
}

double assembleFirePower( const double* attraction ){
    double primary_firepower  [4];
    double secondary_firepower[4];
    double ammo_spend = 0;
    for(const UnitState& unit : composition.units){
        double activity  = unit.n * unit.safeAmmo * unit.supressed;
        double safeAmmo2 = unit.safeAmmo*unit.safeAmmo;
        const double* fps1 = unit.type->primary  ->firepower;
        const double* fps2 = unit.type->secondary->firepower;
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
        double fp2 = primary_firepower[i];
        double vis = attraction[i];
        fp1 *= firePowerSaturation( fp1,     vis );
        fp2 *= firePowerSaturation( fp1+fp2, vis );
        firepower[i] = fp1 + fp2;
    }
    return ammo_spend;
    // choose targets accorgint to attraction
    // if there is not enought visible/attractive targets, we try use secondary weapon
}

void distributeFirePower( const double* firepower){ 
    double normalized_firepower[4];
    for(int i=0; i<4; i++){ 
        normalized_firepower[i] = firepower[i] / target_attraction[i];
        //normalized_firepower[i] = firepower[i] * invAttraction[i]; 
    };
    for(UnitState& unit : composition.units){
        int iclass = unit.type->armorClass;
        double firepower_got = unit.shot_attaction * normalized_firepower[iclass];
        unit.getFire( firepower_got ); // internal damage model
    }
}

};

// ==============================================

struct Combat{
    double dist;
    BattleField conds;
    CombatSide Attacker;
    CombatSide Defender;

void round(){
    Defender.assembleTargetAttraction( conds, dist );
    Defender.assembleTargetAttraction( conds, dist );
    Attacker.assembleFirePower ( Defender.target_attraction );
    Defender.assembleFirePower ( Attacker.target_attraction );
    // distribute firepower and take damage
    Attacker.distributeFirePower( Defender.firepower );
    Defender.distributeFirePower( Attacker.firepower );
    //Attacker.tryAdvance(); // if remain some action points
    //Defender.tryRetreat();
}

};

#endif