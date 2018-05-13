#include "CommandParser.h"
#include "CommandParser2.h"

#include <cstring>

#include "Vec2.h"

class Particle2f{ public:
    Vec2f pos;
    Vec2f vel;

    void move(float dt, const Vec2f& force){
        vel.add_mul( force, dt);
        pos.add_mul( vel,   dt);
        //vel = vel + force * dt;
        //pos = pos + vel   * dt;
    }

    void moveN(int n,float dt){
        Vec2f force;
        float t =0;
        for(int i=0; i<n; i++){
            force   = pos   *-0.5;    // spring
            force.y += -0.5;          // gravity
            force.add_mul(vel,-0.01); // friction
            move(dt, force);
            printf( "%i %f pos=(%f,%f) vel=(%f,%f) \n", i,t, pos.x,pos.y,  vel.x,vel.y  );
            t+=dt;
        };
    };
};

class MyParser : public CommandParser2{ public:
    virtual void* callTable( int ifunc, int nArgs, void** args ){
        void* ret = 0;
        switch(ifunc){
            case 0:{ if(checkFunc(ifunc,"Vec2"))return 0;
                //Vec2f* X = new Vec2f{ *((float*)args[0]), *((float*)args[1]), *((float*)args[2]) };
                Vec2f* X = new Vec2f{ *((double*)args[0]), *((double*)args[1]) };
                //printf( "X=Vec2f(%f,%f)  |  (%f,%f) \n", X->x, X->y,   *((double*)args[0]), *((double*)args[1]) );
                ret=X;
            }break;
            case 1:{ if(checkFunc(ifunc,"Particle"))return 0;
                Particle2f* X = new Particle2f{ *(Vec2f*)args[0], *(Vec2f*)args[1] };
                ret=X;
            }break;
            case 2:{ if(checkFunc(ifunc,"Move"))return 0;
                ((Particle2f*)args[2])->moveN(*(long*)args[0],*(double*)args[1]);
            }break;
        }
        return ret;
    };
};



int main(){

    printf("====== TEST CommandParser1 \n");

    printf(" \n");
    printf("this parser is very simple and robust but is unable to create variables neither evaluate expressions \n");
    printf(" \n");

    CommandParser cmdPars;
    cmdPars.execFile( "data/comands.ini" );

    printf("====== TEST CommandParser2 \n");

    printf(" \n");
    printf("this parser is able to create variables but still not evaluate expressions, there are also no cycles, conditionals and other language structures \n");
    printf(" \n");

    char * code_str = "$pos = Vec2(0.1,2.0); $vel = Vec2(1.0,5.0); $p1 = Particle(pos,vel); Move(5,0.1,p1);";

    MyParser parser;

    parser.registerType    ( "Vec2"     );
    parser.registerType    ( "Particle" );

    parser.registerFunction( "Vec2",     "Vec2",     Strings{"Float", "Float"}           );
    parser.registerFunction( "Particle", "Particle", Strings{"Vec2", "Vec2"}             );
    parser.registerFunction( "Move",     "void",     Strings{"Int", "Float", "Particle"} );

    printf("====== Execute scritp \n");
    parser.execString( 1000, code_str );

    printf("====== Execute Nativ \n");

    Vec2f pos = (Vec2f){0.1,2.0};
    Vec2f vel = (Vec2f){1.0,5.0};
    Particle2f p1  = (Particle2f){pos,vel};
    p1.moveN(5,0.1);


/*
    printf( " %s\n", test_str );
    printf( " ==== \n" );
    parser.parseString( strlen(test_str), test_str );
    //parser.parse( 0, 0, -1 );
    parser.printItemStruct();
*/
}
