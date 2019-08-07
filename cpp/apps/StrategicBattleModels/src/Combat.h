#ifndef Combat_h
#define Combat_h

struct Combat{
    Battlefield* battlefield =0;
    Rival* attacker = 0;
    Rival* defender = 0;
}

void airSuperiority( Air& bf, Rival& att, Rival& def ){

}

void airCombatStep( Battlefield& bf, Rival& att, Rival& def ){
    // TODO:
    // - fighter-fighter combat for air superiority
    //    - Altitude advantage, than speed advantage
    //    - Agility determine probability to be hit by enemy
    // - than fighters try to shoot down bombers
    
    
    
    
}





//
/*

physically based sorties
http://www.wolframalpha.com/input/?i=solve+differential+equation+y%27+%2B+a*y+%2B+b*y%5E2+%3D+0

Force = Gravity + Drag + Thrust

dv(t)/dt = a*v(t) + b*v(t)^2


*/

#endif
