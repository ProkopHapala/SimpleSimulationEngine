
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "SDL_utils.h"
#include "GUI.h"
#include "testUtils.h"

#include "CommodityNetwork.h"

using namespace CommodityNetwork;

int fontTex=0;
char strtemp[1024];

class TestAppCommodityNetwork : public AppSDL2OGL { public:

    std::vector<Vec2d> cityPoss;
    Economy economy;

    bool bRun = false;
    int iStep =0;
    double dt = 0.01;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppCommodityNetwork( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCommodityNetwork::TestAppCommodityNetwork( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    //fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    fontTex     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    GUI_fontTex = fontTex;
    //const nCity = 4;

    //Vec2d* p = new Vec2d{1.0,2.0};

    int nGoods = 8;
    std::string strGoods[nGoods] = {
    //  name &transport_weight, &density, &price_max, &price_normal
        "Lime    1.0 3.3 1000.0   12.0" ,
        "Coal    1.0 1.3 1000.0   33.0" ,
        "IronOre 1.0 2.5 1000.0   95.0",
        "Iron    1.0 7.8 10000.0 300.0" ,
        "Steel   1.0 7.8 10000.0 712.0",

        "Coke    1.0 1.2 2000.0   50.0" ,
        "Pitch   1.0 0.9 10000.0 100.0" ,
        "Sulphur 1.0 2.0 10000.0 150.0" ,

    };
    for(int i=0; i<nGoods; i++){
        /*
        CommodityType* cargo = new CommodityType( strGoods[i].c_str() );
        cargo->id = economy.goods.size();
        economy.goods    .push_back( cargo );
        economy.goodsDict.insert( {cargo->name, cargo} );
        printf( "Comodity %i %s \n", cargo->id, cargo->name.c_str() );
        */
        economy.insert( new CommodityType( strGoods[i].c_str() ) );
    }


    int nTech = 3;
    std::string strTech[nTech] = {
    //             &cycle_time, &unit_space        consumes                  produces       machines
        //"SteelMaking 1.0 1.0 | Iron=10.0  Coal=20.0   |  Steel=9.0   |  people=1"
        //"SteelMaking, 1.0, 1.0| Iron=10.0, Coal=20.0 | Steel=9.0 | people=1"
        //"IronMaking{{1.0;2.0};consume{{IronOre;10.0};{Coal;20.0};{Lime;20.0}};produce{{Iron;9.0}};machines{{people;1}}}",
        //"SteelMaking{{1.0;2.0};consume{{Iron;10.0};{Coal;20.0}};produce{{Steel;9.0}};machines{{people;1}}}"
        "IronMaking{{1.0;2.0};consume{{IronOre;10.0};{Coke;20.0};{Lime;20.0}};produce{{Iron;9.0}};machines{{people;1}}}",
        "SteelMaking{{1.0;2.0};consume{{Iron;10.0};{Coke;20.0}};produce{{Steel;9.0}};machines{{people;1}}}",
        "Coking{{1.0;2.0};consume{{Coal;10.0}};produce{{Coke;7.5};{Pitch;1.0};{Sulphur;0.1}};machines{{people;1}}}"
    };
    for(int i=0; i<nTech; i++){
        /*
        Technology* tech = new Technology( &strTech[i][0], economy.goodsDict );
        tech->id = economy.technologies.size();
        economy.technologies  .push_back( tech );
        economy.technologyDict.insert( {tech->name, tech} );
        printf( "Technology %i %s nConsume %i nPoduce %i \n", tech->id, tech->name.c_str(), tech->consumes.size(), tech->produces.size() );
        */
        economy.insert( new Technology( &strTech[i][0], economy.goodsDict ) );
    }


    //exit(0);
    // make cities
    cityPoss.push_back( {0.0,0.0} );
    cityPoss.push_back( {0.0,1.0} );
    cityPoss.push_back( {1.0,1.0} );
    cityPoss.push_back( {1.0,0.0} );

    const int nRoad = 4;
    Vec2i  links[nRoad] = { {0,1}, {1,2}, {2,3}, {3,0} };

    //for( City* city : economy.cities ){}
    //for(Vec2d& p : cityPoss ){
    //    economy.cities.push_back();
    //};

    // ---- set cities

    DEBUG

    economy.cities.resize( cityPoss.size() );
    for( int ic=0;ic<cityPoss.size();ic++ ){
        City* city = new City();
        economy.cities[ic] = city;
        city->id = ic;
        city->parkCapacity = 1000.0;

        // register each tech in each city
        for( Technology* tech : economy.technologies ){
            printf( "Register Tech %i %s \n", tech->id, tech->name.c_str() );
            city->registerTech(*tech, economy.goods );
        }
        printf( "City %i nGoods %i nTech %i\n", ic, city->goods.size(), city->factories.size() );

        printf("City %i Goods: ", ic ); for(auto it: city->goods){ printf( "{%s %g $%g}", it.second->type->name.c_str(), it.second->ammount, it.second->price ); } printf("\n", ic );

    }

    DEBUG

    // ----- create roads
    economy.roads.reserve(nRoad);
    for(int i=0; i<nRoad; i++){
        Road* rd = new Road();{
            rd->id = i;
            rd->a  = economy.cities[ links[i].a ];
            rd->b  = economy.cities[ links[i].b ];
            rd->length = (cityPoss[links[i].a] - cityPoss[links[i].b]).norm();
        }
        //printf( "  %i %i \n",   links[i].a, links[i].b );
        printf( "road %i %i    %i %i \n", rd->a->id, rd->b->id,    links[i].a, links[i].b );
        economy.roads.push_back( rd );
    }

    DEBUG

     // ----- create cars

    CarType* carType0 = new CarType();{
        carType0->id              = 0;
        carType0->parkSpace       = 1;
        carType0->capacity_weight = 1;
        carType0->capacity_volume = 1;
        carType0->speed           = 1;
    };

    printf( "==== Cars \n" );

    int carPerCity = 1;

    //for( int ic=0;ic<cityPoss.size();ic++ ){
    for( int ird=0;ird<economy.roads.size();ird++ ){
        //City* city = economy.cities[ic] = city;
        Road* rd = economy.roads[ird];
        for(int i=0; i<carPerCity; i++){
            Car* c = new Car();
            c->type = carType0;
            //c->park = city;
            c->road = rd;
            c->tryPark(*rd->a);
            c->fpos = -0.1;
            printf( "Car[%i] road %li (%i,%i)\n", economy.cars.size(),  (long)rd, rd->a->id, rd->b->id );
            economy.cars.push_back(c);
        }
    }

    printf( "==== Cars DONE \n" );

    // ---- set comodity states
    CommodityState* gs;
    //int id_Iron = economy.goodsDict["Iron"]->id;
    //gs = economy.cities[1]->goods[id_Iron];
    gs=economy.getGoodsInCity( 1, "IronOre"); gs->ammount = 200.0;
    gs=economy.getGoodsInCity( 1, "Coal"   ); gs->ammount = 1400.0;
    gs=economy.getGoodsInCity( 1, "Lime"   ); gs->ammount = 600.0;

    //int id_Coal = economy.goodsDict["Coal"]->id;
    //1gs = economy.cities[1]->goods[id_Coal];
    //gs->ammount = 4000.0;
    gs->price   = 30;

    for( int ic=0;ic<economy.cities.size(); ic++ ){
        City* city = economy.cities[ic];
        printf("City %i Goods: ", ic ); for(auto it: city->goods){ printf( "{%s %g $%g}", it.second->type->name.c_str(), it.second->ammount, it.second->price ); } printf("\n", ic );
    }

    DEBUG



    printf( " ===== SETUP DONE ==== \n" );

    DEBUG

} // TestAppCommodityNetwork::TestAppCommodityNetwork




void TestAppCommodityNetwork::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    if(bRun){
        iStep++;
        printf(" ======== iStep %i \n", iStep );
        economy.step( dt );
    }

    for(City* ct: economy.cities){
        //printf( "road %i %i \n", rd->a->id, rd->b->id );
        Draw2D::drawCircle_d(cityPoss[ct->id], 0.1, 16, false );

        sprintf(strtemp, "%i", ct->id );
        //Draw::drawText( strtemp, fontTex, fontSizeDef*0.1, {16,4} );
        Draw2D::drawText( strtemp, cityPoss[ct->id], {2,3}, fontTex, 0.2 );
	}

	for(Road* rd: economy.roads){
        //printf( "road %i %i \n", rd->a->id, rd->b->id );
        Draw2D::drawLine_d( cityPoss[rd->a->id],  cityPoss[rd->b->id] );
	}

    for(Car* c: economy.cars){
        if(c->road){
            Vec2d p; p.set_lincomb( 1-c->fpos, cityPoss[c->road->a->id], c->fpos, cityPoss[c->road->b->id] );
            Draw2D::drawPointCross_d( p, 0.1 );
        }
        //printf( "road %i %i \n", rd->a->id, rd->b->id );
	}

}

void TestAppCommodityNetwork::drawHUD(){
    //char str[1024];
    int ipick = 1;
    economy.cities[ipick]->goodsInfo( strtemp );
    //sprintf( strtemp, "AHOJ !!! \n" );
    //printf( "=== city[%i]\n%s", ct->id, str );
    Draw2D::drawText( strtemp, {100,100}, {160,16*fontSizeDef}, fontTex, fontSizeDef );
}

void TestAppCommodityNetwork::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

void TestAppCommodityNetwork::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppCommodityNetwork::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_SPACE:  bRun = !bRun; break;
            } break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //pickParticle( world.picked );
                    //ipick = pline1.nearestPoint( {mouse_begin_x,mouse_begin_y} );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
                    //ipick = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}



// ===================== MAIN

TestAppCommodityNetwork * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppCommodityNetwork( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















