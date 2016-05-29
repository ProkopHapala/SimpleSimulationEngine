
#ifndef SoldierType_h
#define SoldierType_h

class SoldierType{
	public:
    char   * name = "default";

    // movement
    double   mass      = 1.0;  // collision mass
    double   radius    = 0.5;  // collision radius
    double   max_speed = 1.0;  // speed in good terrain / conditions
    double   min_speed = 1.0;  // speed in bad terrain  / conditions

    // stamina / moral
    double   moral     = 1.0;
    double   stamina   = 1.0;

    //meele combat
    double melee_period      = 1.0;
    double melee_range       = 1.2;
    double melee_damage      = 1.0;
    double melee_penetration = 1.0;
    double melee_dStamina    = 1.0; // stamina cost of meele attack
    double melee_push        = 1.0; // force exerted by weapon in push

    double attack_skill      = 1.0;
    double defence_skill     = 1.0;

    // ranged attack
    double fire_period       = 1.0;
    double fire_range        = 0.0; // if this is <1.0 unit is not ranged;   We may modify this later based on physics
    double fire_damage       = 0.0;
    double fire_penetration  = 0.0;
    double fire_spread       = 0.1; // tg(angular_spread),  spread = distance*fire_spread
    double fire_dStamina     = 1.0; // stamina cost of fire attack

    // armor and defence
    double armorFw           = 1.0;
    double armorBk           = 0.5;

    // shield
    double shield_cos        = 0.5; //
    double shield_miss       = 0.5; // probability to miss the shield
    double shield_endurance  = 1.0;


};

#endif
