
#ifndef CommodityNetwork_h
#define CommodityNetwork_h

// see:  SimpleSimulationEngine/cpp/apps/LandCraft/Economy.h


#include <string>
#include <vector>
#include <unordered_map>

#include "fastmath.h"
#include "str2tree.h"

#include "testUtils.h"



/*

Economy Model 1
----------------
Solve interaction between two basic problems
 * Production - "Factory" is instance of certain "Technology" which converts input{set of goods} to output{set of goods}
 * Transport  - "Car" assigned to certain road transport goods from city "A" to city "B" if cost(B) > cost{A}+transport_cost

*/

static const int nMaxStoreKind = 3;
static int       counter_commodityType = 0;


static double reserveOffset = 100;


namespace CommodityNetwork{

struct Car;
struct City;

struct StoreState{
    int kind=0;
    //double weightMax=0;
    double volume=0;
    double volumeMax=0;
};

struct CommodityType{
    int id=0;
    std::string name;
    int storeKind=0;
    double density         =1;
	double transport_weight=1; // transport weight can difer from unit weight due to container
	double transport_volume=1;
	double price_max       =1;
	double price_normal    =1;

    void fromString(const char* s ){
        char tmp[256]; sscanf( s, "%s %lf %lf %lf %lf", tmp, &transport_weight, &density, &price_max, &price_normal ); name=std::string(tmp);
        transport_volume = transport_weight/density;
    }

    void init0(){
        id=0;
        storeKind=0;
        density         =1;
        transport_weight=1; // transport weight can difer from unit weight due to container
        transport_volume=1;
        price_max       =1;
        price_normal    =1;
    };

    CommodityType() = default;
    CommodityType(const char* s){ fromString(s); storeKind=0; id=counter_commodityType; counter_commodityType++; }
};

struct CommodityState{
    //int id=0;
    CommodityType* type=0;
    double ammount=0;      // current ammount
    double reserve=0;
    double supply =0;      // expected ammount to arrive by import and production
    double demand =0;      // expected ammount to be cosumed
    double price  =1;

    CommodityState( ) = default;
    CommodityState( CommodityType& type_ ){
        type = &type_; ammount=0; reserve=0; supply=0; demand=0; price=type_.price_normal;
    };

    void updatePrice(){
        double surplus = ammount - reserve;
        double change  = supply  - demand;
        double timeFrame = 5.0;
        surplus += change*timeFrame;
        double price_mod = clamp( surplus / (  reserve + reserveOffset ), -1.0, 1.0 );
        price = type->price_normal * ( 1 - 0.5*price_mod );
    }

};

struct Technology{
    int id;
    std::string  name;
    double cycle_time;  // time to produce unit
	double unit_space;  // unit space taken in factory
    std::unordered_map<int,double> consumes;
    std::unordered_map<int,double> produces;
    //std::unordered_map machines;

    /*
    void print(){
        printf("%s %f %f \n", name.c_str(), cycle_time, unit_space );
        //printf("*consumes "); printMap_d(consumes);
        //printf("*produces "); printMap_d(produces);
    }
    */

    void fromString(char* s, std::unordered_map<std::string,CommodityType*>& goodDict ){
        char stmp[256];
        Str2Tree st;
        st.setStr(s);
        printf("%s\n",s);

        //strncpy(s+st.ich,n);
        st.step(stmp); name = stmp; printf("%s\n",stmp);
        st.seekN(1);
        st.step(stmp); cycle_time=atof(stmp); printf("%s\n",stmp);
        st.step(stmp); unit_space=atof(stmp); printf(":%s\n",stmp);
        st.seekN(3);
        int    gid;
        double N;
        do{
            st.step(stmp); gid=goodDict[stmp]->id; printf("<%s   %i \n",stmp,st.level);
            st.step(stmp); N=atof(stmp);           printf(":%s  %i \n",stmp,st.level);
            consumes.insert({gid,N});
            st.seekN(2);
        }while(st.level>1);
        st.seekN(2);
        do{
            st.step(stmp); gid=goodDict[stmp]->id; printf(">%s   %i \n",stmp,st.level);
            st.step(stmp); N=atof(stmp);           printf(":%s  %i \n",stmp,st.level);
            produces.insert({gid,N});
            st.seekN(2);
        }while(st.level>1);
        /*
        st.seekN(2);
        do{
            st.step(stmp); printf(">%s   %i \n",stmp,st.level);
            st.step(stmp); printf(":%s  %i \n",stmp,st.level);
            st.seekN(2);
        }while(st.level>1);
        */
    }

    Technology();
    Technology(char* s, std::unordered_map<std::string,CommodityType*>& goodDict ){fromString(s,goodDict);}

};

struct Factory{
    double space           = 0;
    const Technology* tech = 0; // to do  - why const ?

    /*
    bool setTechnology( Technology* tech_ ){
        //for( auto it : currentTenchnology->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        //for( auto it : currentTenchnology->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        for( auto it : tech->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
        for( auto it : tech->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
        tech = tech;
    }
    */

    double produce( double dt, double Nmax, std::unordered_map<int,CommodityState*>& goods ){
        // check max possible
        /*
        for( auto it : currentTenchnology->consumes ){
            auto got = stored.find(it.first);
            if( got == stored.end() ) return 0;
            //double Nmax_i = got->second->ammount/it.second; //    stored/par_unit
            double Nmax_i = got->second/it.second; //    stored/par_unit
            N=fmin(N,Nmax_i);
        }
        //for( auto it : currentTenchnology->consumes ){ stored[it.first]->ammount -= it.second*N; } // remove resources
        //for( auto it : currentTenchnology->produces ){ stored[it.first]->ammount += it.second*N; } // store products
        for( auto it : currentTenchnology->consumes ){ stored[it.first] -= it.second*N; }
        for( auto it : currentTenchnology->produces ){ stored[it.first] += it.second*N; }
        return N;
        */
        double N = fmin( Nmax, dt*space/tech->unit_space );
        //printf( "produce Nmax %g dt %g space %g/%g N %g \n", Nmax, dt, space, tech->unit_space, N );
        for( auto it : tech->consumes ){ goods[it.first]->ammount -= it.second*N; }
        for( auto it : tech->produces ){ goods[it.first]->ammount += it.second*N; }
        return N; // just place-holder
    }

};

struct TradeStats{
    CommodityState* best_cargo   = 0;
    double          best_profit  = 0;
    double          total_profit = 0;
    void evalRoad( City *A, City *B );
};

struct Road{
    int id;
    City *a=0,*b=0;
    double length=0;
    TradeStats tradeStats[2];

    void evalTradeStats();
    CommodityState* optTransfer(bool bSwap, double maxVol, double& outN, double& best_profit);

};

struct City{
    int id=0;
    //city<>
    double parkUsed    =0;
    double parkCapacity=0;
    double carDemand   =0;
    //std::vector<StoreState>              stores;   // ????? Do We Need this when there is "goods"
    StoreState stores[nMaxStoreKind];
    //std::list<Cars> cars;
    std::unordered_map<int,CommodityState*> goods;     // cargo_id -> store
    std::unordered_map<int,Factory*>        factories; // tech_id  -> Factory
    std::vector<Road*> roads;

    void registerCommodity( CommodityType& type ){
        //printf( "register: %s \n", type.name.c_str() );
        auto got=goods.find(type.id);
        if( got==goods.end() ){
            //CommodityState* s = new CommodityState{ &type, 0, 0, 0, type.price_normal };
            CommodityState* s = new CommodityState(type);
            goods.insert( {type.id, s} );
            //stores[type->storeKind] = ;
            //printf( ".... registered \n" );
            /*
            CommodityType* type=0;
            double ammount=0;      // current ammount
            double supply =0;      // expected ammount to arrive by import and production
            double demand =0;      // expected ammount to be cosumed
            double price  =1;
            */
        }
    }

    /*
    CommodityState* getGoods( const char* name ){
        int id_Iron = goodsDict[name]->id;
        return cities[1]->goods[id_Iron];
    }
    */

    bool registerTech( const Technology& tech, std::vector<CommodityType*>& goodsTypes ){
        auto got = factories.find(tech.id);
        if(got!=factories.end()) return true;
        for( auto it : tech.consumes ){ registerCommodity( *goodsTypes[it.first] ); }
        for( auto it : tech.produces ){ registerCommodity( *goodsTypes[it.first] ); }
        //for( auto it : tech.produces ){ outPrice += goods[it.first]->price; }
        Factory* fac = new Factory();
        fac->space = 10.0;
        fac->tech=&tech;
        factories.insert({tech.id,fac});
        return false;
    }

    double rankTech( const Technology& tech ){
        double inPrice  = 0;
        double outPrice = 0;
        //for( auto it : tech.consumes ){ auto got=goods.find(it.first); if( got!=goods.end() ) inPrice  += got->second->price; }
        //for( auto it : tech.produces ){ auto got=goods.find(it.first); if( got!=goods.end() ) outPrice += got->second->price; }
        // assuming tech is known
        for( auto it : tech.consumes ){ inPrice  += goods[it.first]->price; }
        for( auto it : tech.produces ){ outPrice += goods[it.first]->price; }
        // ToDo : consider storare and factory space
        // ToDo : consider absent commodities
        //printf( "rank price: in(%g) >? out(%g) \n", outPrice, inPrice );
        return outPrice - inPrice; // netGain
    }

    double maxTech( const Technology& tech ){
        double limit = 1e+300;
        for( auto it : tech.consumes ){
            //printf( "consume %s : %g / %g = %g \n", goods[it.first]->type->name.c_str(), goods[it.first]->ammount, it.second, goods[it.first]->ammount/it.second );
            limit = fmin( limit, goods[it.first]->ammount/it.second );
        }
        //printf( "limit N %g \n", limit );
        return limit; // max number units which can be produced using stored input resources
    }

    double produce( double dt ){
        for(auto it: factories){
            Factory* fc = it.second;
            if(fc->space<=0) continue;
            if( rankTech( *fc->tech ) > 0 ){
                double   N = maxTech( *fc->tech );
                fc->produce( dt, N, goods );
            }
        }
    }

    double addCommodity( double N, const CommodityType& cargo ){
        StoreState& ss  = stores[ cargo.storeKind ];
        double dv       = N * cargo.transport_volume;
        //printf( "dv %g \n", dv );
        //printf( "ss.volumeMax %g ss.volume %g \n", ss.volumeMax, ss.volume );
        dv              = fmax( dv, ss.volumeMax - ss.volume ); // ToDo : consider both mass and weight for transport
        //printf( "dv %g \n", dv );
        double da  = dv/cargo.transport_volume;
        ss.volume  +=dv;
        goods[cargo.id]->ammount +=da; // ToDo: check if the commodity is known ?
        return da;
    }

    double takeCommodity( double Nmax, const CommodityType& cargo ){
        //printf( "takeCommodity %s %g %i \n", cargo.name.c_str(), Nmax, cargo.id );
        double N        = fmin( Nmax, goods[cargo.id]->ammount );
        //StoreState& ss  = stores[ cargo.storeKind ];
        //ss.volume      -= N * cargo.transport_volume;
        goods[cargo.id]->ammount -= N; // ToDo: check if the commodity is known ?
        return N;
    }

    int goodsInfo( char* str ){
        char* str0=str;
        str += sprintf( str, "carDemand %g$ \n", carDemand );
        for( auto& it : goods ){
            const CommodityState* gs = it.second;
            str += sprintf( str, "%s$%3.0f %g[kg]\n", gs->type->name.c_str(), gs->price, gs->ammount );
        }
        return str-str0;
    }

    void updatePrices(){
        for( auto& it : goods ){
            it.second->updatePrice();
        }
        double exportProfit = 0;
        for(const Road* rd: roads){
            if(rd->a==this){ exportProfit += rd->tradeStats[0].total_profit; }else{ exportProfit += rd->tradeStats[1].total_profit; }
        }
        carDemand = exportProfit * 0.5;
    }
    //void unload(Car& c){
    //    c.ammount -= addCommodity( c.ammount, *c.cargo );
    //}

};

inline void TradeStats::evalRoad( City *A, City *B ){
    //City *A,*B;
    //if(bSwap){ A=road.b;B=road.a; }else{ A=road.a;B=road.b; }
    best_cargo=0;
    best_profit=0;
    total_profit=0;
    for( auto it : A->goods ){
        auto got = B->goods.find(it.first);
        double Nexport = fmax( 0, it.second->ammount - it.second->reserve );
        if(got!=B->goods.end()){
            double unit_profit = got->second->price - it.second->price;
            if(unit_profit<=0) continue;
            //double maxN   = maxVol/it.second->type->transport_volume;
            //double ammount=0;      // current ammount
            //double reserve=0;
            total_profit  += unit_profit * Nexport;
            if(unit_profit>best_profit){
                //printf( "trade[%i->%i] d$ %g (%g$,%g[kg])->(%g$,%g[kg]) \n", A->id, B->id, unit_profit, it.second->price, it.second->ammount, got->second->price, got->second->ammount );
                best_cargo    = it.second;
                best_profit   = unit_profit;
            }
        }
    }
    //printf( "trade[%i->%i]  d$ %g \n", A->id, B->id, total_profit );
}

inline void Road::evalTradeStats(){
    tradeStats[0].evalRoad( a, b );
    tradeStats[1].evalRoad( b, a );
    //printf( "road[%i][%i->%i] export %g inport %g \n", id, a->id, b->id, tradeStats[0].total_profit, tradeStats[1].total_profit );
};

inline CommodityState* Road::optTransfer(bool bSwap, double maxVol, double& outN, double& best_profit){
    City *A,*B;
    if(bSwap){ A=b;B=a; }else{ A=a;B=b; }
    //int    best_id    =-1;
    CommodityState* best_cargo = 0;
    best_profit= 0;
    //double best_N     = 0;
    //printf( "A,B %i,%i \n", A->id, B->id );
    for( auto it : A->goods ){
        auto got = B->goods.find(it.first);
        if(got!=B->goods.end()){
            //printf( "optTransfer: found %i \n", it.first );
            double profit = got->second->price - it.second->price;
            double maxN   = maxVol/it.second->type->transport_volume;
            double N      = fmin(it.second->ammount, maxN );
            profit *= N;
            //if(profit!=0)printf( "====optTransfer: %s N %g profit %g=%g-%g \n", got->second->type->name.c_str(), N, profit,  got->second->price, it.second->price );
            if(profit>best_profit){
                best_cargo  = it.second;
                best_profit= profit;
                //best_id    = it.first;
                outN       = N;
            }
        }
    }
    //printf( "optTransfer: DONE \n" );
    //printf( "optTransfer: best %li \n", (long)best_cargo );
    //if(best_cargo)printf( "optTransfer: best %s profit %g outN %g \n", best_cargo->type->name.c_str(), best_profit, outN );
    //return id;
    return best_cargo;
}


struct CarType{
    int id;
    double parkSpace       = 1;
    double capacity_weight = 1;
	double capacity_volume = 1;
	double speed           = 1;
	double costPerKm       = 1;
	//double loadTime;
	//double consumption;
};

struct Car{
    int id;
    //int job;
    CarType*              type=0;
    const CommodityType* cargo=0;
    Road*           road=0;
    City*           city=0;
    double       ammount=0;
    double       heading=1;   // direction of transport
    double          fpos=0.5; // position on the road

    virtual void move(double dt){
        //double v = fmax( type->speed, heading );
        double v = type->speed * heading;
        fpos += v*dt/road->length;
    }

    bool tryPark(City& where){
        double space = type->parkSpace;
        if( where.parkUsed + space < where.parkCapacity ){
            //road=0;
            city=&where;
            where.parkUsed+=space;
            return true;
        }
        return false;
    }

    void depart(){
        if      (fpos<0){
            fpos=0; heading=1;
        }else if(fpos>1){
            fpos=1; heading=-1;
        }
    }

    void unload(City* city){
        //printf("car[%i].unload %i\n", id, city->id);
        if(cargo){
            ammount -= city->addCommodity( ammount, *cargo );
            printf( "unload car[%i] city[%i] \n", id, city->id );
        }
        heading=0;
    }
    //void load  (City* city){ ammount += city->takeCommodity( ammount, *cargo ); }

    void tryDepartWithCargo( bool bSwap ){
        double N, profit;
        //printf("loadProfitableCargo swap:%i\n", bSwap );
        double jurneyCost = road->length * type->costPerKm;
        CommodityState* cs = road->optTransfer( bSwap, type->capacity_volume, N, profit );
        //printf( "loadProfitableCargo done road->optTransfer \n" );
        double dCarDemand = road->b->carDemand - road->a->carDemand;
        if(bSwap)dCarDemand=-dCarDemand;
        //printf( "profit %g dCarDemand %g jurneyCost %g \n", profit, dCarDemand, jurneyCost );
        if(dCarDemand>0)profit+=dCarDemand;
        if( profit > jurneyCost ){
            if(cs){
                cargo     = cs->type;
                City* ct  = bSwap?road->b:road->a;
                //printf("loadProfitableCargo city %i \n", ct->id );
                ammount += ct->takeCommodity( N, *cargo );
                //printf("loadProfitableCargo ammount %g \n", ammount );
            }
            depart();
        }else{
            heading=0;
        }
        //printf( "loadProfitableCargo DONE \n" );
    }

    void tryTerminal(){
        //if(heading == 0) return;
        //printf( "Car[%i] City[%i->%i] fpos %g heading %g \n", id, road->a->id, road->b->id, fpos, heading );
        if      ( fpos<=0.0  ){
            if(heading<0){ unload(road->a);           }
            else         { tryDepartWithCargo(false); }
        }else if ( fpos>=1.0 ){
            if(heading>0){ unload(road->b);           }
            else         { tryDepartWithCargo(true);  }
        }
    }

    //void tryChangeRoad( ){}

};


class Economy{ public:

    //std::unordered_map<int,Technology   *> technologies;
    //std::unordered_map<int,CommodityType*> goods;

    std::unordered_map<std::string,Technology   *> technologyDict;
    std::unordered_map<std::string,CommodityType*> goodsDict;

    std::vector<Technology   *> technologies;
    std::vector<CommodityType*> goods;

    std::vector<City*> cities;
    std::vector<Car*>  cars;
    std::vector<Road*> roads;

    int insert( CommodityType* cargo ){
        //CommodityType* cargo = new CommodityType( strGoods[i].c_str() );
        cargo->id = goods.size();
        goods    .push_back( cargo );
        goodsDict.insert( {cargo->name, cargo} );
        printf( "Comodity %i %s \n", cargo->id, cargo->name.c_str() );
        return cargo->id ;
    }

    int insert( Technology* tech ){
        //Technology* tech = new Technology( &strTech[i][0], goodsDict );
        tech->id = technologies.size();
        technologies  .push_back( tech );
        technologyDict.insert( {tech->name, tech} );
        printf( "Technology %i %s nConsume %i nPoduce %i \n", tech->id, tech->name.c_str(), tech->consumes.size(), tech->produces.size() );
        return tech->id;
    }

    void newRoad(int i, int j, double length ){
        Road* rd = new Road();
        rd->id = roads.size();
        rd->a  = cities[ i ];
        rd->b  = cities[ j ];
        roads.push_back(rd);
        cities[ i ]->roads.push_back( rd );
        cities[ j ]->roads.push_back( rd );
        rd->length = length;
        //rd->length = (cityPoss[links[i].a] - cityPoss[links[i].b]).norm();
    }

    void addCar( int iRoad, CarType* carType ){
        Car* c  = new Car();
        c->id   = cars.size();
        c->type = carType;
        //c->park = city;
        Road* rd = roads[iRoad];
        c->road  = rd;
        c->tryPark(*rd->a);
        //c->fpos  = -0.1;
        //printf( "Car[%i] road %li (%i,%i)\n", c->id,  (long)rd, rd->a->id, rd->b->id );
        cars.push_back(c);
    }

    void step(double dt){
        for(City* ct : cities){
            ct->updatePrices();
            ct->produce( dt );
        }
        for( Road* rd: roads){
            rd->evalTradeStats();
        }
        for(Car* c:cars){
            if(c->road){
                c->tryTerminal();
                if( fabs(c->heading)>1e-6 ){
                    c->move(dt);
                }
            }
        }
    }

    CommodityState* getGoodsInCity(int icity, const char* name ){
        return cities[icity]->goods[goodsDict[name]->id];
    }
};


}

#endif

