
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
    //double   moral_regain   = 1.0;
    double   stamina_regain = 0.01;
    double   time_buffer    = 15.0;

    //meele combat
    double melee_period      = 10.0;
    double melee_range       = 1.2;
    double melee_damage      = 3.0;
    double melee_penetration = 1.0;
    double melee_fStamina    = 0.8; // stamina cost of meele attack
    double melee_push        = 1.0; // force exerted by weapon in push

    double attack_skill      = 1.0;
    double defence_skill     = 1.0;
    double defence_period    = 2.0;
    double defence_fStamina  = 0.9;

    // ranged attack
    double fire_period       = 20.0;
    double fire_range        = 0.0; // if this is <1.0 unit is not ranged;   We may modify this later based on physics
    double fire_damage       = 3.0;
    double fire_penetration  = 1.0;
    double fire_spread       = 1.1; // tg(angular_spread),  spread = distance*fire_spread
    double fire_fStamina     = 0.8; // stamina cost of fire attack

    // armor and defence
    double damage_tolerance  = 10.0;
    double armorFw           = 1.0;
    double armorBk           = 0.5;

    // shield
    double shield_cos        = 0.5; //
    double shield_miss       = 0.5; // probability to miss the shield
    double shield_endurance  = 100.0;


};

#endif
