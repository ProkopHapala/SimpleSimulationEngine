
/*

This is combat model for historical line battles ( e.g. Ancient, Medieval, Gunpowder ) ... it is inspired by Paradox games like EU4, CK2, Victoria, Imperator ... 


* units are on discrete rectangular grid
* units can independently decide to advance or retreat
* units preferably attack unit in front of them (same ix), if there no unit they try to attack next nearest unit
* units with exposed flanks tend to panic (some more, some less)


*/

#ifndef SimpleLineBattle_h
#define SimpleLineBattle_h

struct UnitType{

    // bonus matrix ? - bonus against all other units?
    int manuever;
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


}


struct LineCombat{
    int combatWidth;
    Unit* attacker[maxCombatWidth];
    Unit* defender[maxCombatWidth];

    void step(){
        int ixcenter = maxCombatWidth/2;
        for(int ix=0; ix<combatWidth; ix++){
            // --- attacker has initiative
            Unit* att = attacker[ix];
            if( att ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<att->type->manuever; j++ ){
                    Unit* def = defender[ix+j*dx];
                    if( def ){
                        if( def.advance-def.advance ){
                            melee( *att, *def ); 
                            break;
                        } 
                    }
                }
            }
        }
        // --- units which were not engaged by attack
        for(int ix=0; ix<combatWidth; ix++){
            if( def ){
                int dx=(ix<ixcenter)?1:-1;
                for( int j=0; j<att->type->manuever; j++ ){
                    Unit* def = defender[ix+j*dx];
                    if( def ){ melee( *att, *def ); break; }
                }
            }
        }
    }

};





#endif