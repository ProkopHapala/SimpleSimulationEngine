
#ifndef SoldierType_h
#define SoldierType_h

#include <string>
#include <cstdlib>
#include <string.h>
#include <stdio.h>

#include <unordered_map>

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

    char* toStrCaptioned( char * sout ){
        sout += sprintf( sout, "mass              %f\n", mass             );
        sout += sprintf( sout, "radius            %f\n", radius           );
        sout += sprintf( sout, "max_speed         %f\n", min_speed        );
        sout += sprintf( sout, "min_speed         %f\n", min_speed        );
        sout += sprintf( sout, "stamina_regain    %f\n", stamina_regain   );
        sout += sprintf( sout, "time_buffer       %f\n", time_buffer      );
        sout += sprintf( sout, "melee_skill       %f\n", melee_skill      );
        sout += sprintf( sout, "melee_period      %f\n", melee_period     );
        sout += sprintf( sout, "melee_range       %f\n", melee_range      );
        sout += sprintf( sout, "melee_damage      %f\n", melee_damage     );
        sout += sprintf( sout, "melee_penetration %f\n", melee_penetration);
        sout += sprintf( sout, "melee_fStamina    %f\n", melee_fStamina   );
        sout += sprintf( sout, "melee_push        %f\n", melee_push       );
        sout += sprintf( sout, "defence_skill     %f\n", defence_skill    );
        sout += sprintf( sout, "defence_fStamina  %f\n", defence_fStamina );
        sout += sprintf( sout, "fire_period       %f\n", fire_period      );
        sout += sprintf( sout, "fire_range        %f\n", fire_range       );
        sout += sprintf( sout, "fire_damage       %f\n", fire_damage      );
        sout += sprintf( sout, "fire_penetration  %f\n", fire_penetration );
        sout += sprintf( sout, "fire_spread       %f\n", fire_spread      );
        sout += sprintf( sout, "fire_fStamina     %f\n", fire_fStamina    );
        sout += sprintf( sout, "fire_ammo         %i\n", fire_ammo        );
        sout += sprintf( sout, "damage_tolerance  %f\n", damage_tolerance );
        sout += sprintf( sout, "armorFw           %f\n", armorFw          );
        sout += sprintf( sout, "armorBk           %f\n", armorBk          );
        sout += sprintf( sout, "shield_cos        %f\n", shield_cos       );
        sout += sprintf( sout, "shield_miss       %f\n", shield_miss      );
        sout += sprintf( sout, "shield_endurance  %f\n", shield_endurance );
        return sout;
    };

    char* toStrCaptioned_2( char * sout ){
        sout += sprintf( sout, "mass             %f\n", mass             );
        sout += sprintf( sout, "radius           %f\n", radius           );
        sout += sprintf( sout, "max_speed        %f\n", min_speed        );
        sout += sprintf( sout, "min_speed        %f\n", min_speed        );
        sout += sprintf( sout, "stamina_regain   %f\n", stamina_regain   );
        sout += sprintf( sout, "time_buffer      %f\n", time_buffer      );
        sout += sprintf( sout, ":: melee ::\n");
        sout += sprintf( sout, "skill            %f\n", melee_skill      );
        sout += sprintf( sout, "period           %f\n", melee_period     );
        sout += sprintf( sout, "range            %f\n", melee_range      );
        sout += sprintf( sout, "damage           %f\n", melee_damage     );
        sout += sprintf( sout, "penetration      %f\n", melee_penetration);
        sout += sprintf( sout, "fStamina         %f\n", melee_fStamina   );
        sout += sprintf( sout, "push             %f\n", melee_push       );
        sout += sprintf( sout, ":: defence ::\n");
        sout += sprintf( sout, "skill            %f\n", defence_skill    );
        sout += sprintf( sout, "fStamina         %f\n", defence_fStamina );
        sout += sprintf( sout, "demage_tolerance %f\n", damage_tolerance );
        sout += sprintf( sout, "armorFw          %f\n", armorFw          );
        sout += sprintf( sout, "armorBk          %f\n", armorBk          );
        sout += sprintf( sout, "shield_cos       %f\n", shield_cos       );
        sout += sprintf( sout, "shield_miss      %f\n", shield_miss      );
        sout += sprintf( sout, "shield_endurance %f\n", shield_endurance );
        sout += sprintf( sout, ":: fire ::\n");
        sout += sprintf( sout, "period           %f\n", fire_period      );
        sout += sprintf( sout, "range            %f\n", fire_range       );
        sout += sprintf( sout, "damage           %f\n", fire_damage      );
        sout += sprintf( sout, "penetration      %f\n", fire_penetration );
        sout += sprintf( sout, "spread           %f\n", fire_spread      );
        sout += sprintf( sout, "fStamina         %f\n", fire_fStamina    );
        sout += sprintf( sout, "ammo             %i\n", fire_ammo        );
        return sout;
    };

    SoldierType(){};
    SoldierType( char * str ){ fromString( str); };

};

using SoldierTypeDict = std::unordered_map<std::string,SoldierType*>;

#endif
