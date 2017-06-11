
#ifndef SoldierType_h
#define SoldierType_h

#include <string>
#include <cstdlib>
#include <string.h>
#include <stdio.h>

class SoldierType{
	public:
    //char   * name = "default";
    std::string name = "default";

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
    double melee_skill       = 1.0;
    double melee_period      = 10.0;
    double melee_range       = 1.2;
    double melee_damage      = 3.0;
    double melee_penetration = 1.0;
    double melee_fStamina    = 0.8; // stamina cost of meele attack
    double melee_push        = 1.0; // force exerted by weapon in push

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
    int    fire_ammo         = 30;

    // armor and defence
    double damage_tolerance  = 10.0;
    double armorFw           = 1.0;
    double armorBk           = 0.5;

    // shield
    double shield_cos        = 0.5; //
    double shield_miss       = 0.5; // probability to miss the shield
    double shield_endurance  = 100.0;

    //void loadFromFile( char * fname ){};

    void fromString( char * str_ ){
        //printf("1 \n");
        char *str    = strdup(str_);
        printf( "%s\n", str );

        char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );
        name = strdup(token);
        //sscanf( token, "%[^\n]s", name );
        printf( "%s\n", name.c_str() );

        //printf( " basic \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf %lf", &mass, &radius, &max_speed, &min_speed );
        //printf( "%lf %lf %lf %lf\n", mass, radius, max_speed, min_speed );

        //printf( " time \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf", &stamina_regain, &time_buffer );
        //printf( "%lf %lf\n", stamina_regain, time_buffer );

        //printf( " melee \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf %lf %lf %lf %lf", &melee_skill, &melee_period, &melee_range, &melee_damage, &melee_penetration, &melee_fStamina, &melee_push );
        //printf( "%lf %lf %lf %lf %lf %lf %lf %lf\n", melee_skill, melee_period, melee_range, melee_damage, melee_penetration, melee_fStamina, melee_push );

        //printf( " defence \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf", &defence_skill, &defence_period, &defence_fStamina );
        //printf( "%lf %lf %lf %lf\n", defence_skill, defence_period, defence_fStamina );

        //printf( " fire \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf %lf %lf %lf %li", &fire_period, &fire_range, &fire_damage, &fire_penetration, &fire_spread, &fire_fStamina, &fire_ammo );
        //printf( "%lf %lf %lf %lf %lf %lf %li\n", fire_period, fire_range, fire_damage, fire_penetration, fire_spread, fire_fStamina, fire_ammo );

        //printf( " armor \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf", &damage_tolerance, &armorFw, &armorBk );
        //printf( "%lf %lf %lf\n", damage_tolerance, armorFw, armorBk );

        //printf( " shield \n" );
        token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
        sscanf( token, "%lf %lf %lf", &shield_cos, &shield_miss, &shield_endurance );
        //printf( "%lf %lf %lf\n", shield_cos, shield_miss, shield_endurance );

        //exit(0);
    };

    SoldierType(){};
    SoldierType( char * str ){ fromString( str); };

};

#endif
