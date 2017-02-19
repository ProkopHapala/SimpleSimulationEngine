
#include <cstdio>
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

    RigidBody        rb;
    KinematicBody  * kb = &rb;
    PointBody      * pb = &rb;

/*
    //  see IT DOES NOT WORK !!!!
    rb. pos.set(1.0); //printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    //rb._pos::pos.set(1.0); //printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    kb->pos.set(2.0); printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
    pb->pos.set(3.0); printf( "%f %f %f\n", rb.pos.x, kb->pos.x, pb->pos.x );
*/

};
