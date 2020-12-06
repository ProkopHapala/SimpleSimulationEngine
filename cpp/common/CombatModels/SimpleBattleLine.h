
/*

This is combat model for historical line battles ( e.g. Ancient, Medieval, Gunpowder ) ... it is inspired by Paradox games like EU4, CK2, Victoria, Imperator ... 


* units are on discrete rectangular grid
* units can independently decide to advance or retreat
* units preferably attack unit in front of them (same ix), if there no unit they try to attack next nearest unit
* units with exposed flanks tend to panic (some more, some less)


dA/dt = kab * B
dB/dt = kba * A

d_t A/kab = B
d_t B/kba = A

d_t d_t B = kba*kab * B
d_t d_t A = kba*kab * A

kab,kba are damage rates 



Combat Rate -   kba*kab
Combat Bias -   kba/kab ... ratio between killing rate of A, and killing rate of B

k = kba*kab    
f = kba/kab

k*f = kba*kab*kba/kab =  kba*kba = kba^2
kba = sqrt(k*f)
kab = k/kba = kba/f

*/

#ifndef SimpleLineBattle_h
#define SimpleLineBattle_h

struct UnitType{
    double skill;
    double armor;       // [mm]
    double penetration; // [mm]

    // bonus matrix ? - bonus against all other units?
    int manuever;
    int shoot_range;
};

struct Unit{
    UnitType* type;

    bool engaged;
    int advance;
    int n;
    double moral;
};

//const int maxCombatWidth = 5;
const int maxCombatWidth = 40;


void melee( Unit& attacker, Unit& defender ){
    attacker.engaged=true;
    defender.engaged=true;

    double dskill  = ( attacker.type->skill       - defender.type->skill );
    double att_pen = ( attacker.type->penetration - defender.type->armor );
    double def_pen = ( attacker.type->penetration - defender.type->armor );

}

void shoot( Unit& shooter, Unit& target ){


}

struct LineBattle{
    int combatWidth;
    Unit* attacker[maxCombatWidth];
    Unit* defender[maxCombatWidth];

    void step(){
        int ixcenter = maxCombatWidth/2;
        // --- attack
        for(int ix=0; ix<combatWidth; ix++){
            // --- attacker has initiative
            Unit* att = attacker[ix];
            if( att ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<att->type->manuever; j++ ){
                    Unit* def = defender[ix+j*dx];
                    if( def ){
                        if( (def->advance+def->advance)<0 ){
                            melee( *att, *def ); 
                            break;
                        } 
                    }
                }
            }
        }
        // --- units which were not engaged by attack
        for(int ix=0; ix<combatWidth; ix++){
            Unit* def = defender[ix];
            if( def ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<def->type->manuever; j++ ){
                    Unit* att = attacker[ix+j*dx];
                    if( def ){ melee( *def, *att ); break; }
                }
            }
        }
        // shoot phase
        for(int ix=0; ix<combatWidth; ix++){
            Unit* shooter;
            shooter = defender[ix];
            if( shooter ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<shooter->type->shoot_range; j++ ){
                    Unit* target = defender[ix+j*dx];
                    if( target ){ shoot( *shooter, *target ); break; }
                }
            }
            shooter = attacker[ix];
            if( shooter ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<shooter->type->shoot_range; j++ ){
                    Unit* target = defender[ix+j*dx];
                    if( target ){ shoot( *shooter, *target ); break; }
                }
            }
        }
    }

};





#endif