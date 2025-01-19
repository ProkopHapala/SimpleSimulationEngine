
#ifndef  SpaceCraftGUI_h
#define  SpaceCraftGUI_h

#include <unistd.h>
#include <dirent.h>

// #include <SDL2/SDL.h>
// #include <SDL2/SDL_opengl.h>
// #include "Draw.h"
// #include "Draw2D.h"
// #include "Draw3D.h"
// #include "Draw3D_Surf.h"
// #include "SDL_utils.h"

// #include <vector>
// #include "Vec2.h"
// #include "Vec3.h"
// #include "Mat3.h"
// #include "quaternion.h"

//#include "Noise.h"
//#include "SphereSampling.h"
//#include "DrawSphereMap.h"

//#include "SpaceCraft.h"
//#include "EditSpaceCraft.h"
//#include "MeshBuilder.h"

#include "MeshBuilder2.h"
#include "TrussDynamics_f.h"

#include "Interfaces.h"

#include "SpaceCraftComponents.h"
#include "SpaceCraft.h"
#include "GUI.h"

using namespace Mesh;

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

namespace SpaceCrafting{

enum class EDIT_MODE:int{ vertex=1, edge=2, component=4, size }; // http://www.cprogramming.com/c++11/c++11-nullptr-strongly-typed-enum-class.html

class PickerUI{ public:

    EDIT_MODE edit_mode = EDIT_MODE::vertex;
    int picked   =-1;
    //bool bSim = true; // use simulation to pick

    //SpaceCraft*   craft=0;
    //Builder2*  mesh =0;
    //TrussDynamics_f*  sim  =0;

    Picker* picker=0;

    //std::unordered_map<int,GUIPanel*> modePanels; // open proper panel for each mode  ( to simulte behaviour of the switch in SpaceCraftEditorApp::selectCompGui() )
    std::unordered_map<int,BindLoader*> mode_objs;  // open proper panel for each component
    BindLoader* mode_obj;                           // currently selected component 

    double Rmax = 0.0;
    Vec3d ray0;
    Vec3d hray;

    int switch_mode(){
        edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); 
        picked = -1;
        printf("PickerUI::switch_mode() edit_mode=%i\n", (int)edit_mode);
        return (int)edit_mode;
    } 

    void pick(){
        int typ = picker->pick_nearest( ray0, hray, picked, (int)edit_mode, Rmax );
        //printf("picked picked=%i typ=%i edit_mode=%i ray0(%g,%g,%g) hray(%g,%g,%g)\n", picked, typ, (int)edit_mode, ray0.x, ray0.y, ray0.z, hray.x, hray.y, hray.z );
        mode_obj = get( mode_objs, typ, (BindLoader*)0 );
        if(mode_obj) mode_obj->bindLoad( picker->getPickedObject(picked, typ ) );
    }

    void* getPickedObject(){ return picker->getPickedObject(picked, (int)edit_mode ); }

};

class PlateGUI:public BoundGUI{ public:
    //static const int np=4;
    PlateGUI(int xmin, int ymin, int xmax, int dy):BoundGUI("Plate",xmin,ymin,xmax,dy,4){
        opened=false;
        subs[0]->caption="g1.a ";
        subs[1]->caption="g1.b ";
        subs[2]->caption="g2.a ";
        subs[3]->caption="g2.b ";
    }
    int bindLoad(void* o_)override{
        Plate* o = (Plate*) o_;
        //printf( "PlateGUI %li \n", o );
        if(drivers==0)drivers=new GUIPanelWatcher[nsubs];
        drivers[0].bindLoad(subs[0],&o->g1span.a,false);
        drivers[1].bindLoad(subs[1],&o->g1span.b,false);
        drivers[2].bindLoad(subs[2],&o->g2span.a,false);
        drivers[3].bindLoad(subs[3],&o->g2span.b,false);
        binded=true;
        open();
        //redraw=true; tryRender();
        return 0;
    }
};

class GirderGUI:public BoundGUI{ public:
    //static const int np=4;
    GirderGUI(int xmin, int ymin, int xmax, int dy):BoundGUI("Girder",xmin,ymin,xmax,dy,2){
        opened=false;
        subs[0]->caption="nseg ";
        subs[1]->caption="mseg ";
    }
    int bindLoad(void* o_)override{
        Girder* o = (Girder*) o_;
        //printf( "GirderGUI %li \n", o );
        if(drivers==0)drivers=new GUIPanelWatcher[nsubs];
        drivers[0].bindLoad(subs[0],&o->nseg,true);
        drivers[1].bindLoad(subs[1],&o->mseg,true);
        binded=true;
        open();
        //redraw=true; tryRender();
        return 0;
    }
};

} // namespace SpaceCrafting

#endif
