
#ifndef CommodityNetwork_h
#define CommodityNetwork_h

// see:  SimpleSimulationEngine/cpp/apps/LandCraft/Economy.h


#include <string>
#include <vector>
#include <unordered_map>

#include "fastmath.h"
#include "str2tree.h"


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
    CommodityType();
    CommodityType(const char* s){fromString(s);}
};

struct CommodityState{
    //int id=0;
    CommodityType* type=0;
    double ammount=0;      // current ammount
    double supply =0;      // expected ammount to arrive by import and production
    double demand =0;      // expected ammount to be cosumed
    double price  =1;

    CommodityState( );
    CommodityState( CommodityType& type_ ){
        type = &type_; ammount=0; supply=0; demand=0; price=type_.price_normal;
    };
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
    const Technology* tech = 0;

    /*
    bool setTechnology( Technology* tech_ ){
        //for( auto it : currentTenchnology->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        //for( auto it : currentTenchnology->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        for( auto it : tech->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
        for( auto it : tech->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
        tech = tech;
    }
    */

    double produce(double N){
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
        return -1; // just place-holder
    }

};


struct City{
    int id=0;
    //city<>
    double parkUsed    =0;
    double parkCapacity=0;
    std::vector<StoreState>                 stores;
    //std::list<Cars> cars;
    std::unordered_map<int,CommodityState*> goods;     // cargo_id -> store
    std::unordered_map<int,Factory*>        factories; // tech_id  -> Factory

    void registerCommodity( CommodityType& type ){
        auto got=goods.find(type.id);
        if( got==goods.end() ){
            //CommodityState* s = new CommodityState{ &type, 0, 0, 0, type.price_normal };
            CommodityState* s = new CommodityState(type);
            goods.insert( {type.id, s} );
            /*
            CommodityType* type=0;
            double ammount=0;      // current ammount
            double supply =0;      // expected ammount to arrive by import and production
            double demand =0;      // expected ammount to be cosumed
            double price  =1;
            */
        }
    }

    bool registerTech( const Technology& tech, std::vector<CommodityType*>& goodsTypes ){
        auto got = factories.find(tech.id);
        if(got!=factories.end()) return true;
        for( auto it : tech.consumes ){ registerCommodity( *goodsTypes[it.first] ); }
        for( auto it : tech.produces ){ registerCommodity( *goodsTypes[it.first] ); }
        //for( auto it : tech.produces ){ outPrice += goods[it.first]->price; }
        Factory* fac = new Factory();
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
        return outPrice - inPrice;
    }

    double addCommodity( double N, const CommodityType& cargo ){
        StoreState& ss  = stores[ cargo.storeKind ];
        double dv       = N * cargo.transport_volume;
        dv              = fmax( dv, ss.volumeMax - ss.volume );
        double da  = dv/cargo.transport_volume;
        ss.volume  +=dv;

        goods[cargo.id]->ammount +=da; // ToDo: check if the commodity is known ?

        /*
        auto got = goods.find( cargo.id );
        if ( got != goods.end() ){
            got->second->ammount +=da;
        }else{
            // ToDo - is this even possible ?
            //        import or produce unknown commodity ?
            CommodityState* snew = new CommodityState();
            {
                //snew.id  =
                snew->type    = &cargo;
                snew->ammount = da;
                snew->demand  = 0;
                snew->price   = cargo.price_normal;
            }
            goods.insert( {cargo.id,snew} );
        }
        */
        return da;
    }

    //void unload(Car& c){
    //    c.ammount -= addCommodity( c.ammount, *c.cargo );
    //}

};

struct Road{
    int id;
    City *a=0,*b=0;
    double length=0;

    void rankTransfer(){
        for( auto it : a->goods ){
            auto got = b->goods.find(it.first);
            if(got!=b->goods.end()){

            }
        }
    }


};

struct CarType{
    int id;
    double parkSpace       = 1;
    double capacity_weight = 1;
	double capacity_volume = 1;
	double speed           = 1;
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
    double       heading=0;   // direction of transport
    double          fpos=0.5; // position on the road

    virtual void move(double dt){
        double v = fmax( type->speed, heading );
        fpos += v*dt/road->length;
    }

    bool tryPark(City& where){
        double space = type->parkSpace;
        if( where.parkUsed + space < where.parkCapacity ){
            road=0;
            city=&where;
            where.parkUsed+=space;
            return true;
        }
        return false;
    }

    void unload(City* city){ ammount -= city->addCommodity( ammount, *cargo ); }

    void tryUndload(){
        if (heading > 0){
            if ( fpos>1.0 ){
                unload(road->b);
            }
        }else{
            if( fpos<0.0 ){
                unload(road->b);
            }
        }
    }

    void tryLoad(){

    }

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

    void step(double dt){
        for(Car* c:cars){
            if(c->city) c->tryLoad();
            if(c->road) c->move(dt);
            if(c->city) c->tryUndload();
        }
        //for(Car* c:cars){
        //    c->tryUndload();
        //}
    }

};


}

#endif

