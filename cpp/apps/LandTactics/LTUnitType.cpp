#include "LTUnitType.h" // THE HEADER



//////////////////////////////////////////
////            LTGunType
//////////////////////////////////////////

void LTGunType::updateAux(){
    crossArea   = caliber*caliber*0.25*M_PI;
    balisicCoef = crossArea/pMass;
    //AP = ; // ?
};

void LTGunType::fromString(const char * str_ ){
    //printf("1 \n");
    char *str    = strdup(str_);
    //printf( "%s\n", str );
    char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );
    name = strdup(stripWhite(token));
    //sscanf( token, "%[^\n]s", name );
    printf( "name: >>%s<<\n", name.c_str() );
    //printf( " basic \n" );
    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf %li %lf %lf", &rps, &nburst, &spread, &vMuzzle, &pMass );
    printf( "rps, nburst, spread, vMuzzle, pMass %lf %li %lf %lf %lf\n", rps, nburst, spread, vMuzzle, pMass );

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf %lf", &AP, &ExDMG );
    printf( "AP, ExDMG %lf %lf\n", AP, ExDMG );

    token = strtok(NULL, ";");
    dmgType=str2enum( token, 4, sDmgTypes );
    printf( " dmgType: >>%s<< %i \n", token, dmgType );

};

char* LTGunType::toStrCaptioned( char * sout ){
    sout += sprintf( sout, "== GunType    %s\n", name.c_str()  );
    sout += sprintf( sout, "Nburst        %i\n", nburst        );
    sout += sprintf( sout, "round/min.    %f\n", nburst*rps*60 );
    sout += sprintf( sout, "acc.(1/tg)    %f\n", 1/spread      );
    sout += sprintf( sout, "vMuzzle[m/s]  %f\n", vMuzzle       );
    sout += sprintf( sout, "pMass[kg]     %f\n", pMass         );
    sout += sprintf( sout, "AP[mm]        %f\n", AP*1000       );
    sout += sprintf( sout, "KE.Dmg. [kJ]  %f\n", getKineticDamage(0,0) );
    sout += sprintf( sout, "Ext.Dmg.[kJ]  %f\n", ExDMG         );
    sout += sprintf( sout, "DmgType       %s\n", sDmgTypes[dmgType] );
    return sout;
};


//////////////////////////////////////////
////            LTUnitType
//////////////////////////////////////////

void  LTUnitType::updateAux(){
    szAreas.x = sz.y*sz.z;
    szAreas.y = sz.x*sz.z;
    szAreas.z = sz.x*sz.y;
};

void LTUnitType::fromString(const char * str_, GunTypeDict& gunTypeDict ){
    //printf("1 \n");
    char *str    = strdup(str_);
    //printf( "%s\n", str );
    char * token = strtok(str, ";"); //printf( ">>%s<<\n", token );
    name = strdup(stripWhite(token));
    //sscanf( token, "%[^\n]s", name );
    printf( "name: >>%s<<\n", name.c_str() );

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    kind  = str2enum( token, 5, sUnitKind );
    printf( " kind: >>%s<< %i \n", token, kind );

    //printf( " basic \n" );
    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf %lf %lf",  &sz.a, &sz.b, &sz.c );
    printf( "size:  %lf %lf %lf\n", sz.a,  sz.b,  sz.c );
    //printf( " time \n" );
    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token,    "%lf %lf %lf",  &mass, &maxSpeed, &enginePower );
    printf( "mobility: %lf %lf %lf\n", mass,  maxSpeed,  enginePower);
    //printf( " melee \n" );
    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf %lf %lf %lf %lf",  &armorFront, &armorSide, &armorBack, &armorTop, &armorBottom );
    printf( "armor: %lf %lf %lf %lf %lf\n", armorFront,  armorSide,  armorBack,  armorTop,  armorBottom );
    //printf( " defence \n" );

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%lf",  &HP );
    printf(        "%lf\n", HP );

    token = strtok(NULL, ";"); //printf( ">>%s<<\n", token );
    sscanf( token, "%li\n", &nGun );

    for(int i=0; i<nGun; i++){
        token = strtok(NULL, ";");
        std::string s = stripWhite( token );
        auto it = gunTypeDict.find( s );
        if(it!=gunTypeDict.end()){
            printf("found >>%s<<\n", s.c_str() );
            guns[i] = it->second;
            printf( "gun[%i].name >>%s<< \n", i, guns[i]->name.c_str() );
        }else{
            printf("not found >>%s<<\n", s.c_str() );
            exit(-1);
        };
    }
    updateAux();
};

char* LTUnitType::toStrCaptioned( char * sout, bool bGunDetials ){
    sout += sprintf( sout, "=== UnitType     %s\n", name.c_str() );
    sout += sprintf( sout, "size       [m]   %f %f %f\n", sz.a,sz.b,sz.c );
    sout += sprintf( sout, "mass       [kg]  %f\n", mass              );
    sout += sprintf( sout, "maxSpeed   [m/s] %f\n", maxSpeed          );
    sout += sprintf( sout, "enginePower[kW]  %f\n", enginePower*0.001 );
    sout += sprintf( sout, "HitPoints        %f\n", HP                );
    sout += sprintf( sout, ":: Armor ::\n"  );
    sout += sprintf( sout, "front   %f\n", armorFront  );
    sout += sprintf( sout, "side    %f\n", armorSide   );
    sout += sprintf( sout, "back    %f\n", armorBack   );
    sout += sprintf( sout, "top     %f\n", armorTop    );
    sout += sprintf( sout, "bottom  %f\n", armorBottom );
    sout += sprintf( sout, "Side    %f\n", armorSide   );
    sout += sprintf( sout, "Guns :: \n" );
    //printf( "DEBUG 1 \n" );
    for(int i=0; i<nGun; i++){
        sout += sprintf( sout, "%s, ", guns[i]->name.c_str() );
    }
    sout += sprintf( sout, "\n" );
    if( bGunDetials ){
        for(int i=0; i<nGun; i++){
            //printf( "DEBUG ::%i \n", i );
            //sout += sprintf( sout, "gun[%i] %s : \n", i, guns[i]->name.c_str() );
            sout  = guns[i]->toStrCaptioned( sout );
        }
    }
    return sout;
};
