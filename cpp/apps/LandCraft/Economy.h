#ifndef Economy_h
#define Economy_h

#include <string>
#include <vector>
#include <unordered_map>

void printMap_d( const std::unordered_map<std::string,double>& mapa ){
    for(auto it : mapa ){ printf("%s=%f ", it.first.c_str(), it.second ); }; printf("\n");
}

void str2map( char * str, std::unordered_map<std::string,double>& mapa ){
    //char s[256];
    //strcpy(s, "one two three");
    //printf( "%s\n", str );
    char* tok = strtok(str, " \t");
    while (true){
        //printf("token: %s\n", token);
        tok = strtok(NULL, " \t");
        if(tok==NULL) break;
        char* sep = strchr( tok,'=');
        if(sep==NULL) break;
        sep[0]='\0';
        double f;
        sscanf(sep+1,"%lf", &f );
        //printf( " %s = %f \n", std::string(tok).c_str(), f );
        mapa[std::string(tok)] = f;
        //printf( " %s = %f \n", tok, f );
    }
    //printMap_d( mapa );
}

class Commodity{ public:
    std::string name;
	double transport_weight; // transport weight can difer from unit weight due to container
	double transport_volume;
	double price_max;
	double price_normal;
    void fromString( char* s ){ char tmp[256]; sscanf( s, "%s %lf %lf %lf %lf", tmp, &transport_weight, &transport_volume, &price_max, &price_normal ); name=std::string(tmp);}
};

class Technology{ public:
    std::string   name;
    double cycle_time;  // time to produce unit
	double unit_space;  // unit space taken in factory
    std::unordered_map<std::string,double> consumes;
    std::unordered_map<std::string,double> produces;
    //std::unordered_map machines;

    bool fromFile( FILE* pFile){
        char tmp1[4096];
        char tmp2[256];
        char *str;
        char c;
        do{
            str = fgets( tmp1, 4096, pFile );
            if(str==NULL) return false;
            c=str[0];
        }while( !(((c>='a')&&(c<='z'))||((c>='A')&&(c<='Z'))) );
        sscanf( str, "%s %lf %lf %lf %lf", tmp2, &cycle_time, &unit_space ); name=std::string(tmp2);
        //printf("%s %f %f \n", name.c_str(), cycle_time, unit_space );
        //exit(0);
        str = fgets( tmp1, 4096, pFile ); str2map( str, consumes );
        str = fgets( tmp1, 4096, pFile ); str2map( str, produces );
        return true;
    }

    void print(){
        printf("%s %f %f \n", name.c_str(), cycle_time, unit_space );
        printf("*consumes "); printMap_d(consumes);
        printf("*produces "); printMap_d(produces);
    }

};

/*
class Store{ public:
    Commodity* commodity = NULL;
    double     ammount   = NULL;
    Store(){};
    Store(Commodity* commodity_,double ammount_){ commodity=commodity_; ammount=ammount_; };
};
*/

class Factory{ public:
    //std::unordered_map<std::string,Store*> stored;
    std::unordered_map<std::string,double> stored;
    Technology* currentTenchnology = NULL;

    void setTechnology( Technology* tech ){
        currentTenchnology = tech;
        //for( auto it : currentTenchnology->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        //for( auto it : currentTenchnology->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert( {it.first, new Store(NULL,0.0)} ); }
        for( auto it : currentTenchnology->consumes ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
        for( auto it : currentTenchnology->produces ){ auto got=stored.find(it.first); if( got==stored.end() ) stored.insert({it.first, 0.0}); }
    }

    double produce(double N){
        // check max possible
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
    }
};

/*
class Factory{ public:
    technology
};
*/

#endif
