
#ifndef spaceCraftEditorUtils_h
#define spaceCraftEditorUtils_h

#include "Tree.h"
#include "Truss.h"

//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

#include <unistd.h>
#include <dirent.h>


int dir2tree(TreeViewTree& node, char * name, const std::string& prefix="", bool bPrint=false ){
    node.content.caption = name;
    std::string path;
    if (prefix.length()==0){ path = name;              } 
    else                   { path= (prefix+"/")+name;  }
    DIR           *dir=0;
    struct dirent *ent=0;
    //if( chdir(name)==0 ){
    if( (dir = opendir( path.c_str() ))!=0){
        if(bPrint)printf("dir2tree(%s)\n", path.c_str() );
        while( (ent = readdir(dir)) != 0){
            //printf("dir '%s' \n", path.c_str() );
            if((ent->d_name[0]=='.'))continue;
            TreeViewTree* tr = new TreeViewTree();
            tr->parrent = &node;
            node.branches.push_back( tr );
            dir2tree( *node.branches.back(), ent->d_name, path );
        }
        closedir(dir);
    }else{
        if(bPrint)printf("leaf '%s'\n", path.c_str() );
    }
    return 0;
}


int makeTruss( Truss& truss ){

    //truss.girder1( (Vec3d){-5.0,0.0,0.0}, (Vec3d){5.0,0.0,0.0}, (Vec3d){0.0,1.0,0.0}, 5, 1.0 );

    Truss trussPlan;
    trussPlan.loadXYZ(  "data/octShip.xyz" );
    //trussPlan.affineTransform( (Mat3d){5.5,0.0,0.0, 0.0,5.5,0.0, 0.0,0.0,5.5}, false );
    trussPlan.affineTransform( (Mat3d){6,0.0,0.0, 0.0,6,0.0, 0.0,0.0,6}, false );
    GirderParams * gpar = new GirderParams [trussPlan.edges.size()];
    //Vec3d ups  = new Vec3d[trussPlan.edges.size()];
    Vec3d ups[] = {
        (Vec3d){0.0,-1.0,0.0},
        (Vec3d){0.0,+1.0,0.0},
        (Vec3d){0.0,0.0,-1.0},
        (Vec3d){0.0,0.0,+1.0},
        (Vec3d){-1.0,0.0,0.0},
        (Vec3d){+1.0,0.0,0.0}
    };
    //truss.makeGriders( 6, &trussPlan.edges[0], &trussPlan.points[0], gpar, ups );
    //truss.makeGriders( trussPlan, gpar, ups, NULL );
    truss.wheel( {0.0,0.0,0.0}, {10.0,0.0,0.0}, {0.0,1.0,0.0}, 50, 0.5 );
    std::vector<Vec2i> ends;
    truss.makeGriders( trussPlan, gpar, ups, &ends );
    truss.autoBridge(ends.size(), &ends[0], 0.8, 0 );
    truss.panel( {0.0,0.0,0.0}, {5.0,0.0,0.0}, {0.0,5.0,0.0}, {5.0,5.0,0.0}, {5,5}, 1.0 );
    delete [] gpar;

    int glo = glGenLists(1);
    glNewList( glo, GL_COMPILE );
    SpaceCrafting::drawTruss( truss, true );
    glEndList();

    return glo;

}

#endif
