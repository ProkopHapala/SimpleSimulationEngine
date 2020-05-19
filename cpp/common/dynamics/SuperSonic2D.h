
#ifndef SuperSonic2D_h
#define SuperSonic2D_h

/*

https://en.wikipedia.org/wiki/Oblique_shock
http://www.i-asem.org/publication_conf/anbre15/T2I.5.AS502_1358F1.pdf
http://www.tfd.chalmers.se/~nian/courses/compflow/notes/TME085_L05.pdf

*  tangential velocity component does not change across the shock
It is always possible to convert an oblique shock into a normal shock by a Galilean transformation.

tg(theta) = 2*cotg(beta) ( (M^2)*sin(beta)^2 -1 ) / ( M^2 (kapa + cos(2*beta) ) +2 )

sin(theta)/cos(theta) =    2*(cos(beta)/sin(beta))   ( (M^2)*sin(beta)^2 -1 ) / ( M^2 (kapa + cos(2*beta) ) +2 )

st/ct  = 2* (cb/sb)  =  ( (M^2)*sb^2 -1 ) / ( M^2 (k + cb^2 - sb^2) ) +2 )


NASA Oblique Shock calculator
https://www.grc.nasa.gov/WWW/K-12/airplane/oblique.html
https://web.archive.org/web/20121021100737/http://www.aerostudents.com/files/aerodynamicsC/obliqueShockWaves.pdf

https://www.youtube.com/watch?v=JvrWmg8m5is


*/



#include "fastmath.h"
#include "Vec2.h"
//#include "Vec3.h"



struct obliqueShock{
    double Min;    // incident mach number
    Vec2d  normal; // direction of shock


};




struct SScell{
    Vec2d  v; // velocity x,y
    double p; // pressure
};

class SuperSonic2D{ public:

    Vec2i ns;
    double** lines;
    double*  line;
    Vec2d* v;


    void eval(){
        // push neighboring samples asside to fullfill equation of state and continuity equation
        //pV = nRT
        for(int i=0; i<ns.y;  i++){
        }
    }



};

#endif

