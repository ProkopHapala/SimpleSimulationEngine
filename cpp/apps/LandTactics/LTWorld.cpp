
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "LTFaction.h"
#include "LTWorld.h" // THE HEADER

#define DEBUG_PLOT_INTERACTION( pa, pb, R, G, B ) if(interacts){ \
    glColor3f( R, G, B ); \
    Draw2D::drawLine_d( pa->pos, pb->pos ); };

void LTWorld::update( ){
    for( int i=0; i<per_frame; i++  ){
        simulationStep( dt );
    }
}

void LTWorld::simulationStep( double dt ){
    for( LTSquad* f : squads ){
        f->update( dt );
    }
};


int LTWorld::getUnitAt( const Vec2d& p, LTFaction * faction ){
    int imin=-1;
    double r2min = 1e+300;
    //for( Unit * u : squads ){ // this does not seem to work
    for( int i=0; i<squads.size(); i++ ){
        LTSquad * u = squads[i];
        if( u->faction == faction ) continue;
        double r2 = p.dist2( u->pos );
        if( (r2<sq(u->radius))&&(r2 < r2min) ){ r2min=r2; imin=i; }
    }
    //printf( " imin %i r2min %f \n", imin, r2min );
    return imin;
};

void LTWorld::initStaticObject(){

    int n      = 128;
    int maxTry = 16;
    //objects = new LTStaticObject[n];
    Rect2d span = (Rect2d){ (Vec2d){map_center.x-200.0,map_center.y-200.0}, (Vec2d){map_center.x+200.0,map_center.y+200.0} };

    square_ruler.setup( {0.0,0.0}, {100.0,100.0} );
    square_ruler.setN ( {64,64} );
    squares = new LTMapSquare[ square_ruler.ntot ];

    int otiles[4];
    for(int i=0; i<n; i++){
        //objects.push_back( LTStaticObject() );
        //LTStaticObject o = objects.back();
        LTStaticObject o;

        o.type = &objectTypes[0]; // TODO: change in future
        //o.id     = i;
        //o->kind   = LTSObjKind::tree;
        o.radius = 10.0;
        //double ri2 = o.radius*o.radius;
        bool pass = false;
        for( int itry=0; itry<maxTry; itry++ ){  // check overlap
            o.pos.set( randf(span.x0,span.x1), randf(span.y0,span.y1) );
            pass = true;
            for( LTStaticObject& oj : objects ){ // TODO: this is brute-force ... use square_ruler to accelerate it
                Vec2d d; d.set_sub( oj.pos, o.pos );
                double r2 = d.norm2();
                double Rmin = oj.radius+o.radius;
                if( r2 < (Rmin*Rmin) ){ pass=false; break; };
            }
            //printf( "   itry %i pass: %i \n", itry, pass);
            if(pass) break;
        }
        //printf( "obj %i pass: %i \n", i, pass);
        if(!pass) break;

        float phi = randf(0.0,M_PI*2);
        o.dir.set( cos(phi), sin(phi) );

        objects.push_back( o );

        int nret = square_ruler.getOverlapingTiles( o.pos, o.radius, otiles );
        for( int j=0; j<nret; j++ ){ squares[ otiles[j] ].objects.push_back(&o); }
    }

    printf( "LTWorld::initStaticObject DONE \n" );
}

int LTWorld::loadGunTypes(const char * fname ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    while ( fgets( str , nbuff, pFile) != NULL ){
        printf("===str: >>%s<<\n", str);
        if (str[0]=='#') continue;
        LTGunType t = LTGunType(str);
        gunTypes.push_back( t );
    }
    fclose(pFile);
    return gunTypes.size();
}

int LTWorld::loadUnitTypes(const char * fname ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    while ( fgets ( str , nbuff, pFile) != NULL ){
        printf("==>>%s<<\n", str);
        if (str[0]=='#') continue;
        LTUnitType t = LTUnitType(str, gunTypeDict);
        unitTypes.push_back( t );
    }
    fclose(pFile);
    return unitTypes.size();
}

/*
int LTWorld::loadSquads(const char * fname ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    while ( fgets ( str , nbuff, pFile) != NULL ){
        printf("==>>%s<<\n", str);
        if (str[0]=='#') continue;
        LTUnitType t = LTUnitType(str, gunTypeDict);
        unitTypes.push_back( t );
    }
    fclose(pFile);
    return unitTypes.size();
}
*/

void LTWorld::init(){
    printf( "========== LTWorld::init() \n" );
    evalAuxSimParams();

    ruler.setSize(128,128);
    ruler.setStep(50);

    map_center = (Vec2d){ruler.na*0.75*ruler.step,ruler.nb*0.5*ruler.step};

    ground = new double[ruler.ntot];
    hydraulics.setSize(ruler.na,ruler.nb);
    hydraulics.ground = ground;

    //world.hydraulics.allocate( 512, 512 );

    hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );
    for( int j=0; j<500; j++ ){
        int isz = 25;
        int ix0 = rand()%(hydraulics.nx-isz);
        int iy0 = rand()%(hydraulics.ny-isz);
        hydraulics.errodeDroples( 200, 100, 0.02, 0.15, 0.5, ix0, iy0, ix0+isz, iy0+isz );
    }
    for(int i=0; i<ruler.ntot; i++){ ground[i] *= maxHeight; };

    //for(int i=0; i<ruler.ntot; i++){ ground[i] = randf(0.0,500.0); };
    //for(int ib=0; ib<ruler.nb; ib++){  for(int ia=0; ia<ruler.na; ia++){  ground[ib*ruler.na+ia] = ia/(float)ruler.na;  } };

    //soldierTypes.push_back( SoldierType(){"pikemen",1.0d,0.25d,1.0d} );
    //soldierTypes.push_back( {"pikemen",1.0d,0.25d,1.0d, 1.0, 1.0 } );
    //unitTypes.push_back( UnitType() );

    // init Object Types




    //loadGunTypes ("data/GunTypes.ini");
    //load2vector( "data/GunTypes.ini", gunTypes );
    processFileLines( "data/GunTypes.ini", [&](char* s){ LTGunType t(s); gunTypes.push_back(t); } );
    vec2map( gunTypes, gunTypeDict );
    printDict( gunTypeDict  );
    //exit(0);
    printf("DEBUG ==== loadGunTypes DONE \n");
    //loadUnitTypes("data/UnitTypes.ini");
    //load2vectorDict( "data/UnitTypes.ini", unitTypes, gunTypeDict );
    processFileLines( "data/UnitTypes.ini", [&](char* s){ LTUnitType t(s,gunTypeDict); unitTypes.push_back(t); } );

    vec2map( unitTypes, unitTypeDict );
    printf("DEBUG ==== loadUnitTypes DONE \n");


    //LTObjectType* ot = new LTObjectType();
    objectTypes.push_back(LTRectHouseType());
    //LTObjectType* ot = new LTRectHouseType();
    LTObjectType* ot = &objectTypes[0]; // TODO change in future
    ot->render_glo();
    //objectTypes.push_back(ot);

    initStaticObject();

    LTSquad * u;

    Vec2d rot,pLook;
    rot = (Vec2d){0.0,+1.0};
    pLook = map_center+(Vec2d){0.0,300.0};

    LTFaction* fac1 = new LTFaction( "RedArmy" , 0xFF0080ff );   factions.push_back( fac1 );
    LTFaction* fac2 = new LTFaction( "BlueArmy", 0xFFff8000 );   factions.push_back( fac2 );

    LTFaction* fac;
    auto lamb_squadLine = [&](char* s) {
        LTSquad* u=new LTSquad(s,unitTypeDict);
        fac->squads.push_back(u);
        squads.push_back(u);
        u->faction=fac;
        u->pos.add(map_center);
        u->populate(u->n);
    };
    fac=fac1; processFileLines( "data/Faction1_squads.ini", lamb_squadLine );
    fac=fac2; processFileLines( "data/Faction2_squads.ini", lamb_squadLine );

    //LTUnitType* ut = unitTypes; // TODO : initialize units from file

    //load2vectorDict( "data/Faction1_squads.ini", fac1->squads, unitTypeDict );
    //load2vectorDict( "data/Faction2_squads.ini", fac2->squads, unitTypeDict );

    /*
    u = new LTSquad( &unitTypes[0], fac1, map_center+(Vec2d){-50.0,-30.0} ); fac1->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;
    u = new LTSquad( &unitTypes[0], fac1, map_center+(Vec2d){ 0.0, -30.0} ); fac1->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;
    u = new LTSquad( &unitTypes[0], fac1, map_center+(Vec2d){+50.0,-30.0} ); fac1->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;
    printf("DEBUG 3 \n");
    rot = (Vec2d){0.0,-1.0};
    pLook = map_center+(Vec2d){0.0,-300.0};

    u = new LTSquad( &unitTypes[0], fac2, map_center+(Vec2d){-50.0, 30.0} ); fac2->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;
    u = new LTSquad( &unitTypes[0], fac2, map_center+(Vec2d){ 00.0, 30.0} ); fac2->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;
    u = new LTSquad( &unitTypes[0], fac2, map_center+(Vec2d){+50.0, 30.0} ); fac2->squads.push_back(u); squads.push_back(u);  u->lookAt(pLook); u->rot=rot;

    for( LTSquad* s : squads ){
        s->populate( 5 );
    }
    */

    printf( "========== LTWorld::init() DONE \n" );
};






