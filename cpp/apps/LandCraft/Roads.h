#ifndef Roads_h
#define Roads_h

#include <string>
//#include <vector>
#include <list>
//#include <unordered_map>
#include "Economy.h"

class RoadTile{ public:
    uint16_t ia,ib;
    double height;
    void print(){printf( "(%i,%i) %f\n", ia, ib, height );};
    //int state; e.g. damaged
};

class Road{ public:
    int  n=0;
    //int* path = NULL;
    RoadTile* path = NULL;
    //std::vector<RoadTile> path;
    //std::vector<int>   path;
    //std::vector<RoadVehicle*> vehicles;
};

class RoadBuilder{ public:    // makes road easy to edit
    //std::list<int> path;
    std::list<RoadTile> path;
    Road* road;

    void pushStright( const Vec2i& ip ){
        RoadTile& back = path.back();
        back.print();
        uint16_t ia=back.ia;
        uint16_t ib=back.ib;
        int sia=1,sib=1;
        //printf( "=========== pushStright (%i,%i) (%i,%i)\n", ia, ib, ip.a, ip.b );
        int dia = ia - ip.a; if( dia<0 ){ int sia = -1; dia=-dia; }
        int dib = ib - ip.b; if( dib<0 ){ int sia = -1; dib=-dib; }
        for(int i=0; i<dia; i++){
            ia += sia;
            path.push_back( {ia,ib,0.0} );
        }
        for(int i=0; i<dib; i++){
            ib += sib;
            path.push_back( {ia,ib,0.0} );
        }
    }

    void readIt(){
        path.clear();
        for(int i=0; i<road->n; i++){
            path.push_back( road->path[i] );
        }
    }

    void writeIt(){
        //Road* road = new Road();
        if( road->n != path.size() ){
            if( road->path ) delete road->path;
            road->n = path.size();
            road->path = new RoadTile[ road->n ];
        }
        int i=0;
        for (auto const& tile : path ){
            road->path[i] = tile;
            i++;
        }
    }
};

class RoadVehicleType{ public:
    double maxSpeed = 50.0;      // [m/s]
    double maxForce = 10.0e+3;   // [N]
    double power    = 100.0e+3;  // [W]
    double mass     = 1e+3;      // [kg]
    double cargoMass   = 1.0e+3; // [kg]
    double cargoVolume = 5.0;    // [m^3]

    double getSpeed( double slope ){ return 1; };
};

class RoadVehicle{ public:
    Road            * road = NULL;
    RoadVehicleType * type = NULL;
    //int  itile=0;
    int  ipath=0;
    //bool dir=1;
    int idir=1;
    double t_rest=0.0;
    bool onWay = false;

    void depart(){
        if(ipath==0      ) idir= 1;
        if(ipath==road->n) idir=-1;
        onWay=true;
    }

    bool moveStep(double& t_left ){
        t_left += t_rest;
        //std::vector<RoadTile>& path = road->path;
        RoadTile* path = road->path;
        //double slope = heights[itile+idir] - heights[itile]; // maybe use road->heights[ipath]
        double slope = path[ipath+idir].height - path[ipath].height;
        double tStep = type->getSpeed( slope );
        double dt = t_left - tStep;
        //printf( "%i %i %f %f\n", ipath, idir, t_left, t_rest );
        if(dt>0){
            if( idir==1 ){
                ipath++; if( ipath==road->n ) return true;
            }else{
                ipath--; if( ipath==0       ) return true;
            }
            t_left -= tStep;
            t_rest = 0;
        }else{
            t_left = 0;
            t_rest = t_left;
        }
        return false;
    }

    void move(double t){
        if( onWay ){
            while( t>0 ){
                if( moveStep( t ) ){
                    //printf( "arrived \n" );
                    onWay = false;
                    break;
                };
            }
        }
    }
};

#endif
