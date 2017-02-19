

/*
with virtual inheritance
_pos      : 3.382 ticks/call ( 3.38172e+08 1e+08 ) | 5.00007e+07
PointBody : 13.658 ticks/call ( 1.36577e+09 1e+08 ) | 5.00007e+07
RigidBody : 44.158 ticks/call ( 4.41584e+09 1e+08 ) | 5.00007e+07


class PointBody : public  _pos, public Updateable {

_pos      : 3.435 ticks/call ( 3.43501e+08 1e+08 ) | 5.00007e+07
PointBody : 22.626 ticks/call ( 2.26263e+09 1e+08 ) | 5.00007e+07
RigidBody : 35.059 ticks/call ( 3.50588e+09 1e+08 ) | 5.00007e+07

*/


#include <cstdio>
#include "testUtils.h"

#include "Bodies.h"

KinematicBody   a0;
PointBody       a1;
RigidBody       a2;
PassiveObject   a3;
ActiveObject    a4;
AttachedObject  a5;

int main(){
    int sz_Collidable = sizeof(Collidable); printf(  "Collidable  %i [Byte]\n", sz_Collidable );
    int sz_Updateable = sizeof(Updateable); printf(  "Updateable  %i [Byte]\n", sz_Updateable );

    int sz__identifiers = sizeof(_identifiers); printf(  "_identifiers  %i [Byte]\n", sz__identifiers );
    int sz__pos         = sizeof(_pos );        printf(  "_pos          %i [Byte]\n", sz__pos );
    int sz__rot         = sizeof(_rot);         printf(  "_rot          %i [Byte]\n", sz__rot );
    int sz__size        = sizeof(_size);         printf(  "_size         %i [Byte]\n", sz__size  );
    int sz__attached    = sizeof(_attached);    printf(  "_attached     %i [Byte]\n", sz__attached );

    int sz_KinematicBody  = sizeof(KinematicBody);  printf(  "KinematicBody   %i [Byte]\n", sz_KinematicBody  );
    int sz_PointBody      = sizeof(PointBody );     printf(  "PointBody       %i [Byte]\n", sz_PointBody  );
    int sz_RigidBody      = sizeof(RigidBody);      printf(  "RigidBody       %i [Byte]\n", sz_RigidBody  );
    int sz_PassiveObject  = sizeof(PassiveObject);  printf(  "PassiveObject   %i [Byte]\n", sz_PassiveObject  );
    int sz_ActiveObject   = sizeof(ActiveObject);   printf(  "ActiveObject    %i [Byte]\n", sz_ActiveObject  );
    int sz_AttachedObject = sizeof(AttachedObject); printf(  "AttachedObject  %i [Byte]\n", sz_AttachedObject  );


    printf(  "RigidBody - componenets %i [Byte]\n", sz_RigidBody - sz_PointBody - sz__rot - sizeof(Mat3d) - 2*sizeof(Vec3d)  );

    /*
    RigidBody        rb;
    KinematicBody  * kb = &rb;
    PointBody      * pb = &rb;


    //  see IT DOES NOT WORK !!!!
    rb. pos.set(1.0); printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    //rb._pos::pos.set(1.0); //printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    kb->pos.set(2.0); printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    pb->pos.set(3.0); printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    */


    constexpr int n = 1000000;
    _pos      * pos_ = new _pos[n];
    PointBody * pbs  = new PointBody[n];
    RigidBody * rbs  = new RigidBody[n];
    RigidBody_flat * rbfs = new RigidBody_flat[n];


    for(int i=0; i<n; i++){
        double x = randf();
        pos_[i].pos.x = x;
        pbs [i].pos.x = x;
        rbs [i].pos.x = x;
        rbfs [i].pos.x = x;
    }

    int m = 100;
    int ncall = n*m;
    double sum;
    long tstart,time;

    sum = 0.0;
    tstart = getCPUticks();
    for(int j=0; j<m; j++){ for(int i=0; i<n; i++){ sum  += pos_[i].pos.x; } }
    time   = getCPUticks() - tstart;
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", "_pos     ", time/double(ncall), (double)time, (double)ncall, sum );

    sum = 0.0;
    tstart = getCPUticks();
    for(int j=0; j<m; j++){ for(int i=0; i<n; i++){ sum  += pbs[i].pos.x; } }
    time   = getCPUticks() - tstart;
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", "PointBody", time/double(ncall), (double)time, (double)ncall, sum );

    sum = 0.0;
    tstart = getCPUticks();
    for(int j=0; j<m; j++){ for(int i=0; i<n; i++){ sum  += rbs[i].pos.x; } }
    time   = getCPUticks() - tstart;
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", "RigidBody", time/double(ncall), (double)time, (double)ncall, sum );

    sum = 0.0;
    tstart = getCPUticks();
    for(int j=0; j<m; j++){ for(int i=0; i<n; i++){ sum  += rbfs[i].pos.x; } }
    time   = getCPUticks() - tstart;
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", "RigidBody_flat", time/double(ncall), (double)time, (double)ncall, sum );

};
