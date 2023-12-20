
#include "Warrior25D.h" // THE HEADER

void Warrior25D::fromString( char const* str ){
    char meshname[256];
    sscanf( str, "%s %lf %lf %lf %lf\n", meshname, &span.c, &mass, &I, &power );
    printf(      "%s %lf %lf %lf %lf\n", meshname, span.c,  mass,  I,  power );
    setMass( mass );
    setI   ( I    );
    OMesh * mesh_ = new OMesh();
    int res = mesh_->fromFileOBJ( meshname );
    if(res>0) mesh=mesh_;
}
