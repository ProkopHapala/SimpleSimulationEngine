
#ifndef MechGrid2D_h
#define MechGrid2D_h

#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"

class MechGrid2D{ public:

/*
    Multi-material mechanical ractangular grid
    * neighboruing pixel with same material are treated as usual (fluid simulation)
    * Boundary between materials is treated by s_edge parameter which define triangle within volume
 *--------*
 |     /  |
 |    /   |
 |   /    |
 |  /     |
 *--------*
*/

    double edge_lenght = 1.0;

    Vec2i ns;

    // Vert properties
    int*    imat     =0;

    // Volume properties
    double* pressure =0;
    double* mass     =0;
    Vec2d*  momentum =0;

    // edge properties
    Vec2d* s_edge = 0; // position of boundary between materials on the edge
    Vec2d* force;



    void pressure2force(){
        int i=0;
        for(int iy=1;iy<ns.y; iy++){
            for(int ix=1;i<ns.x; ix++){
                force[i].x += (pressure[i] - pressure[i+1   ])*edge_lenght;
                force[i].y += (pressure[i] - pressure[i+ns.x])*edge_lenght;
                i++;
            }
        }
    }

    void updateMomnetum(double dt){
        int ntot = ns.x*ns.y;
        for(int i=0;i<ntot; i++){
            Vec2d f = { force[i].x + force[i+1].x, force[i].y + force[i+ns.x].y };
            momentum[i].add_mul( f, dt );
            // need to update edge boundaries
        }
    }

    void advec(double dt){
        int ntot = ns.x*ns.y;
        for(int i=0;i<ntot; i++){
            //s_edge[]
        }
    }

};

#endif

