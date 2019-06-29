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
    ~Road(){ if(path) delete [] path; }
};

class RoadBuilder{ public:    // makes road easy to edit
    //std::list<int> path;
    std::list<RoadTile> path;
    Road* road;

    void pushStright( const Vec2i& ip ){
        RoadTile& back = path.back();
        /*
        back.print();
        uint16_t ia=back.ia;
        uint16_t ib=back.ib;
        int sia=1,sib=1;
        //printf( "=========== pushStright (%i,%i) (%i,%i)\n", ia, ib, ip.a, ip.b );
        //printf( "p0(%i,%i) p2(%i,%i) d(%i,%i) \n", ia,ib,   ip.a,ip.b );
        int dia = ip.a - ia; if( dia<0 ){ sia = -1; dia=-dia; }
        int dib = ip.b - ib; if( dib<0 ){ sib = -1; dib=-dib; }

        if( (sia<0)&&(sib<0) ){ sia };

        //printf( "s(%i,%i) d(%i,%i) \n", sia,sib, dia, dib );
        for(int i=0; i<dia; i++){
            ia += sia;
            path.push_back( {ia,ib,0.0} );
        }
        for(int i=0; i<dib; i++){
            ib += sib;
            path.push_back( {ia,ib,0.0} );
        }
        */

        int dx = ip.x - back.ia;
        int dy = ip.y - back.ib;
        // From http://www-cs-students.stanford.edu/~amitp/Articles/HexLOS.html
        // assume abs(dx) >= abs(dy), it's symmetric otherwise
        int sig,factor,x,y,xone,yone;
        // this is (2); the next line changes from "==" to "!=" if
        // hexagons are not stacked vertically (see first paragraph)
        //sig = ( sign(dx) == sign(dy) );
        if(dx<0){ xone = -1;} else{ xone = 1; } // unit increments
        if(dy<0){ yone = -1;} else{ yone = 1; } // unit increments
        sig = ( xone != yone );



        if (dx % 2)  {  // so dx is divisible by two
            dx *= 2;
            dy *= 2;
        }
        int xerr=dx,yerr=dy;
        dx = abs(dx);
        dy = abs(dy); // don't need the signs anymore

        const bool xbranch = (dx>=dy);
        if( xbranch ){ factor = dx/2; }else{  factor = dy/2; };
        // should start at 0.5, which is just (dx/2)/dx
        x = back.ia;
        y = back.ib;
        //process(x,y);
        path.push_back( {x,y,0.0} );
        int iter=0;
        //while (x != ip.x || y != ip.y) {
        int err=xerr*xerr+yerr*yerr;
        while ( xerr || yerr ) {
            if(xbranch){
                factor += dy;
                if (factor >= dx) {
                    // Make a "diagonal move" in grid (ie vertical or horizontal)
                    factor -= dx;
                    if(sig) {  // vertical move: 2 moves in 1
                        x += xone; // add 1 in the appropriate sign
                        y += yone;
                    } else {   // horizontal move: 2 moves in 2
                        x += xone; // doesn't matter which increases first
                        //process(x,y);
                        path.push_back( {x,y,0.0} );
                        y += yone;
                    }
                } else {
                    x += xone;
                }
            }else{ // y_branch
                factor += dx;
                if (factor >= dy) {
                    // Make a "diagonal move" in grid (ie vertical or horizontal)
                    factor -= dy;
                    if(sig) {  // vertical move: 2 moves in 1
                        x += xone; // add 1 in the appropriate sign
                        y += yone;
                    } else {   // horizontal move: 2 moves in 2
                        y += yone; // doesn't matter which increases first
                        //process(x,y);
                        path.push_back( {x,y,0.0} );
                        x += xone;
                    }
                } else {
                    y += yone;
                }
            };
            //process(x,y);
            path.push_back( {x,y,0.0} );
            printf("iter %i x,y %i,%i    -> %i,%i \n", iter, x,y, ip.x, ip.y   );
            xerr = x - ip.x;
            yerr = y - ip.y;
            int err_=xerr*xerr+yerr*yerr;
            if( (iter>100)||(err_>err) ) break;
            err=err_;
            iter++;
        }

        //exit(0);
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
