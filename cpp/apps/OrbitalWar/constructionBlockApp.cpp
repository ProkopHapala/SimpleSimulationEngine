#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "globals.h"

// Math & Geometry
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Solids.h"
#include "raytrace.h"

// Graphics & GUI
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"

#include "MeshBuilder2.h"
#include "DrawUV.h"
#include "ConstructionBlock.h"
#include "MeshBuilder2Draw.h"
#include "argparse.h"


//#include "testUtils.h"

// ======================  Global Variables & Declarations

Mesh::Builder2 truss;      // final truss to be rendered
Mesh::Builder2 skelet;    // low resolution mesh to start from
ConstructionBlock block;
BlockBuilder bb;    

const int GIRDER_TYPE = 1;
const int ROPE_TYPE   = 2;

Vec3d pivot_point{ 5.0, 0.0, 3.0 };

using namespace Mesh;

void createSkelet(Mesh::Builder2& mesh, std::vector<double>& node_sizes) {
    mesh.clear();
    node_sizes.clear();
    // Nodes from ship_ICF_marksman_2.lua
    node_sizes.push_back(1.0);  mesh.vert({0.0,  0.0, 0.0 }); // 0
    node_sizes.push_back(0.5);  mesh.vert({0.0,  0.0,-15.0}); // 1
    node_sizes.push_back(0.5);  mesh.vert({0.0,  0.0, 20.0}); // 2
    node_sizes.push_back(0.25); mesh.vert({0.0,-10.0, 0.0 }); // 3
    node_sizes.push_back(0.25); mesh.vert({0.0, 10.0, 0.0 }); // 4
    node_sizes.push_back(0.25); mesh.vert({-10.0,0.0, 0.0 }); // 5
    node_sizes.push_back(0.25); mesh.vert({10.0, 0.0, 0.0 }); // 6
    // Girders
    mesh.edge(0, 1, GIRDER_TYPE, -1);
    mesh.edge(0, 2, GIRDER_TYPE, -1);
    mesh.edge(0, 3, GIRDER_TYPE, -1);
    mesh.edge(0, 4, GIRDER_TYPE, -1);
    mesh.edge(0, 5, GIRDER_TYPE, -1);
    mesh.edge(0, 6, GIRDER_TYPE, -1);
    // Ropes
    mesh.edge(1, 3, ROPE_TYPE, -1); mesh.edge(1, 4, ROPE_TYPE, -1); mesh.edge(1, 5, ROPE_TYPE, -1); mesh.edge(1, 6, ROPE_TYPE, -1);
    mesh.edge(2, 3, ROPE_TYPE, -1); mesh.edge(2, 4, ROPE_TYPE, -1); mesh.edge(2, 5, ROPE_TYPE, -1); mesh.edge(2, 6, ROPE_TYPE, -1);
    mesh.edge(3, 5, ROPE_TYPE, -1); mesh.edge(5, 4, ROPE_TYPE, -1); mesh.edge(4, 6, ROPE_TYPE, -1); mesh.edge(6, 3, ROPE_TYPE, -1);
}

void trussFromSkelet(const Mesh::Builder2& skelet, const std::vector<double>& node_sizes, BlockBuilder& bb, Mesh::Builder2& truss) {
    // skelet (Mesh::Builder2) -> BlockBuilder
    bb.clear();
    for (int i = 0; i < skelet.verts.size(); ++i) { bb.addBlock(skelet.verts[i].pos, node_sizes[i]);}
    for (const auto& edge : skelet.edges)         { if (edge.w == GIRDER_TYPE) { bb.connectBlocks(edge.x, edge.y); }}
    Mesh::ConstructionBlockToMeshBuilder cbm;
    // BlockBuilder -> truss (Mesh::Builder2)
    cbm.mesh = &truss;
    cbm.drawBlockBuilder(bb, 4);
    // ropes
    for (const auto& edge : skelet.edges) {
        if (edge.w == ROPE_TYPE) { truss.rope(skelet.verts[edge.x].pos, skelet.verts[edge.y].pos, -1); } // -1 for auto-segmentation
    }
}

void testConstructionBlocks(BlockBuilder& bb, Mesh::Builder2& mesh){
    Mesh::ConstructionBlockToMeshBuilder cbm;
    cbm.mesh = &mesh;
    bb.clear();
    int ic  = bb.addBlock( Vec3d{0.0,  0.0, 0.0 }, 1.0  );
    int ibk = bb.addBlock( Vec3d{0.0,  0.0,-15.0}, 0.5  );
    int ifw = bb.addBlock( Vec3d{0.0,  0.0, 20.0}, 0.5  );
    int ilf = bb.addBlock( Vec3d{0.0,-10.0, 0.0 }, 0.25 );
    int irt = bb.addBlock( Vec3d{0.0, 10.0, 0.0 }, 0.25 );
    bb.connectBlocks(ic,ibk);
    bb.connectBlocks(ic,ifw);
    bb.connectBlocks(ic,ilf);
    bb.connectBlocks(ic,irt);
    cbm.drawBlockBuilder( bb, 4 );
}






// ====================== Class Definitions

class ConstructionBlockApp : public AppSDL2OGL_3D { public:

    int perFrame = 10;

	bool bSmoothLines = 1;
	bool bWireframe   = 1;

    // View control properties
    bool bViewSkelet       = false;
    bool bViewBlockBuilder = false;
    bool bViewMesh         = true;
    bool bViewEdges        = true;
    bool bViewTris         = true;
    bool bViewFaces        = true;
    bool bViewFaceNormals  = false;
    bool bViewPointLabels  = false;
    bool bViewFaceLabels   = false;
    bool bViewTriLabels    = false;
    bool bViewPivotPoint   = true;
    bool bViewSelection    = true;

	//DropDownList lstLuaFiles;
    GUI gui;
    CheckBoxList* viewControls = 0;
    MultiPanel*   contextMenu  = 0;

    //EDIT_MODE edit_mode = EDIT_MODE::component;
    //int picked = -1;
    //Vec3d hray; //= (Vec3d)(cam.rot.c);
    //Vec3d ray0; //= (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);

    int picked_block = -1;
    int ipick = -1;

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( );
	//virtual void camera();
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys ) override;
    //virtual void mouseHandling( );

    void initGUI();

	ConstructionBlockApp( int& id, int WIDTH_, int HEIGHT_ , int argc, char *argv[]);

};

ConstructionBlockApp::ConstructionBlockApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    //Lua1.init();
    fontTex       = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    GUI_fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    Draw::fontTex = fontTex;

    // --- Create truss by extruding Parabola  (e.g. to build magnetic nozzle for nuclear spacecraft)
    //Parabola2Mesh(mesh,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  true );
    //Parabola2Mesh(mesh,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  false ); // does not work - crash in Mesh::drawFace
    //Parabola_ExtrudedWire( mesh, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.5, 0.1 );

    // --- Create truss from low resolution skeleton
    // std::vector<double> node_sizes;
    // createSkelet(skelet, node_sizes);
    // trussFromSkelet(skelet, node_sizes, bb, truss);
    // printf("skelet.printSizes(): "); skelet.printSizes();
    // printf("truss.printSizes():  "); truss.printSizes();
    // skelet.write_obj("skelet.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
    // truss .write_obj("truss.obj",  (uint8_t)(ObjMask::Verts | ObjMask::Edges | ObjMask::Tris));
    
    // --- Create truss using construction block manually
    // testConstructionBlocks( bb, truss);
    // printf("truss.printSizes():         \n"); truss.printSizes();
    // truss.write_obj("high_res_truss.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges | ObjMask::Tris));

    // --- Single block testing
    block.Ls=Vec3d{1.1,1.0,0.9};
    for(int i=0;i<6;i++){
        block.faces[i].typ=1;
    }

    truss.printSizes();
    

    //printf( "truss.tris.size(): \n", truss.tris.size() );

    //plateGui  = (PlateGUI* )gui.addPanel( new PlateGUI ( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );
    //girderGui = (GirderGUI*)gui.addPanel( new GirderGUI( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );

    //truss.loadXYZ(  "data/octShip.xyz" );
    //DropDownList* lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    //lstLuaFiles = new DropDownList( "lua files",20,HEIGHT_-100,200,5); gui.addPanel(lstLuaFiles);
    //lstLuaFiles->setCommand([this](GUIAbstractPanel* panel){ onSelectLuaShipScript.GUIcallback(panel); });

    truss.bAdditiveSelect = true;
    truss.selection_mode = (int)Mesh::Builder2::SelectionMode::vert;
    //mesh.selection_mode = (int)Mesh::Builder2::SelectionMode::edge; break;
    //mesh.selection_mode = (int)Mesh::Builder2::SelectionMode::face; break;
    initGUI();
}

void ConstructionBlockApp::initGUI(){
    // Setup GUI checkboxes for view control
    viewControls = new CheckBoxList();
    viewControls->caption = "View Controls";
    viewControls->initCheckBoxList(5, 5, 150);
    viewControls->addBox("Low-Res Graph", &bViewSkelet);
    viewControls->addBox("Block Builder", &bViewBlockBuilder);
    viewControls->addBox("Mesh",          &bViewMesh);
    viewControls->addBox("Edges",         &bViewEdges);
    viewControls->addBox("Face Normals",  &bViewFaceNormals);
    viewControls->addBox("Point Labels",  &bViewPointLabels);
    viewControls->addBox("Face Labels",   &bViewFaceLabels);
    viewControls->addBox("Tri Labels",    &bViewTriLabels);
    viewControls->addBox("Pivot Point",   &bViewPivotPoint);
    gui.addPanel(viewControls);

    //mp= new MultiPanel( "Edit", gx.x0, ylay.x0, gx.x1, 0,-13); 

    // context menu for right mouse button
    contextMenu = new MultiPanel( "Context Menu", 0, 0, fontSizeDef*20, fontSizeDef*2, -1 );
    contextMenu->hideOnCommand = true;
    //contextMenu->addPanel( "select along line",  {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->ffl.print_nonbonded();   return 0; }; 
    contextMenu->addButton( "selectVertsAlongPolyline", [&](GUIAbstractPanel* p){ 
        printf( "selectVertsAlongPolyline \n" );
        printf( "select along line BEFORE \n" ); truss.printSelectedVerts();
        truss.selectVertsAlongPolyline( 0.1, true ); 
        printf( "select along line AFTER \n" ); truss.printSelectedVerts();
        return 0; 
    } );
    contextMenu->addButton( "plateBetweenEdges", [&](GUIAbstractPanel* p){ 
        printf( "plateBetweenEdges \n" );
        printf( "plateBetweenEdges BEFORE \n" ); truss.printSelectedVerts();
        truss.plateBetweenEdges();
        printf( "plateBetweenEdges AFTER \n" ); truss.printSelectedVerts();
        return 0; 
    } );
    gui.addPanel( contextMenu );
}

void ConstructionBlockApp::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.8f, 0.8f, 0.8f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	//glDisable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
    glLightfv( GL_LIGHT0, GL_POSITION,  (float*)&cam.rot.c  );

    if(bViewSkelet) {
        glDisable(GL_LIGHTING);
        glColor3f(0.8, 0.2, 0.2);
        glLineWidth(3.0);
        drawEdges(skelet);
        glPointSize(10.0);
        drawVerts(skelet);
        glEnable(GL_LIGHTING);
    }

    if(bViewBlockBuilder) {
        Draw3D::drawBlockBuilder( bb );
    }

    if(bViewMesh) {
        glEnable(GL_LIGHTING);
        drawMesh      ( truss, bViewFaces, bViewTris, bViewFaceNormals, bViewEdges);
        drawMeshLabels( truss, bViewFaces, bViewPointLabels, bViewFaceLabels, bViewTriLabels);
        glLineWidth(1.0);
        glColor3f(0.0,0.7,0.0);
        if(bDragging && iDraggingButton==SDL_BUTTON_LEFT ){ drawMuseSelectionBox(); }
        glDisable(GL_DEPTH_TEST);
        if(bViewSelection) {
            drawMeshSelection(truss, ipick);
        }
    }
    if(bViewPivotPoint) { Draw3D::drawPointCross( pivot_point, 0.5 ); }
    //glLineWidth(5.0); Draw3D::drawAxis(10.0);
};

void ConstructionBlockApp::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);

    // void Draw::drawText( const char * str, int itex, float sz, Vec2i block_size ){

    // sprintf(str_tmp, "time=%10.5f[s] mass=%g cog(%g,%g,%g) vcog(%g,%g,%g) L(%g,%g,%g) torq(%g,%g,%g) |F|=%g \n", sim.time, sim.mass, sim.cog.x,sim.cog.y,sim.cog.z, sim.vcog.x,sim.vcog.y,sim.vcog.z, sim.L.x,sim.L.y,sim.L.z, sim.torq.x,sim.torq.y,sim.torq.z, sim.F_residual );
    // //sprintf( str_tmp, "time= %10.5f[s] \n ", sim.time );
    // Draw::drawText( str_tmp, fontTex, fontSizeDef,  {WIDTH,HEIGHT-20}  );

    gui.draw();
    //glColor3f(1.0f,1.0f,1.0f);   txtStatic.view3D( {5,5}, fontTex, 8 );
    //glPopMatrix();
}

void ConstructionBlockApp::keyStateHandling( const Uint8 *keys ){
    //Mat3d camMat;
    //qCamera.toMatrix_T(camMat);
    //Draw3D::drawMatInPos(  (Mat3f)camMat, (Vec3d){0.0,0.0,0.0} );
    //qCamera.toMatrix(camMat);
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.roll(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.roll( -keyRotSpeed ); }

    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b, +0.05*zoom ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -0.05*zoom );  }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a, +0.05*zoom );  }
    //if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ theSpaceCraft->nodes[7]->calong-=0.001; theSpaceCraft->nodes[7]->updateBound(); }
    //if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ theSpaceCraft->nodes[7]->calong+=0.001; theSpaceCraft->nodes[7]->updateBound(); }



};

void ConstructionBlockApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    if(event.type == SDL_MOUSEWHEEL){
        if     (event.wheel.y > 0){ zoom/=VIEW_ZOOM_STEP; }
        else if(event.wheel.y < 0){ zoom*=VIEW_ZOOM_STEP; }
    }
    gui.onEvent(mouseX,mouseY,event);
    switch( event.type ){

         // ========  Mouse events  ========

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    //printf("SDL_BUTTON_LEFT: mouse_ray0 %g %g %g mouse_begin %g %g \n", ray0.x, ray0.y, ray0.z, mouse_begin_x, mouse_begin_y); 
                    mouseStartDragging( SDL_BUTTON_LEFT ); 
                    //picker.pick();
                } break;
                case SDL_BUTTON_RIGHT:{
                    mouseStartDragging( SDL_BUTTON_RIGHT ); 
                } break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( ray0.dist2(ray0_start)<0.1 ){ // too small for selection box 
                        ipick = truss.pickSelect( (Vec3d)ray0, (Vec3d)cam.rot.c, 0.1 );
                        printf( "ipick %i \n", ipick );
                        //if(mesh.bAdditiveSelect){ mesh.selection.push_back(ipick); mesh.printSelection(); }
                    }else{
                        //ipick=-1;
                        int nsel = truss.selectRect( (Vec3d)ray0_start, (Vec3d)ray0, (Mat3d)cam.rot );
                        if(nsel==0) truss.clearSelection();
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT: { 
                    if( ray0.dist2(ray0_start)<0.1 ){ // right click, not dragging
                        contextMenu->showAsContextMenu( mouseX, mouseY );
                        ipick=-1;
                    } 
                } break;
            }
            break;

        // ========  Keyboard events ========

        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_l:{} break;
                case SDLK_f:{
                    truss.selectionToFace();
                } break;

                case SDLK_e:{
                    if( (ipick>=0) && ( truss.selection_mode == (int)Mesh::Builder2::SelectionMode::face ) ){
                        truss.extrudeFace( ipick, 5.0 );    
                    }
                } break;

                // ---- Selection mode  I, O, P
                case SDLK_i: truss.selection_mode = (int)Mesh::Builder2::SelectionMode::edge; truss.bAdditiveSelect = true;  break;
                case SDLK_o: truss.selection_mode = (int)Mesh::Builder2::SelectionMode::vert; truss.bAdditiveSelect = true;  break;
                case SDLK_p: truss.selection_mode = (int)Mesh::Builder2::SelectionMode::face; truss.bAdditiveSelect = false; break;

            }
            break;
    };
    AppSDL2OGL::eventHandling( event );

}

// ===================== MAIN

int main(int argc, char *argv[]){    
    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	ConstructionBlockApp * app = new ConstructionBlockApp( junk , dm.w-150, dm.h-100, argc, argv );
    //app->bindSimulators( &W ); 

    LambdaDict funcs;
    //funcs["-s"]={1,[&](const char** ss){ reloadShip( ss[0] ); }}; 
    //funcs["-perframe"]={1,[&](const char** ss){            sscanf( ss[0], "%i", &app->perFrame );              printf( "COMMAND LINE: -perframe(%i) \n", app->perFrame ); } };
    //funcs["-method"  ]={1,[&](const char** ss){ int im;    sscanf( ss[0], "%i", &im );  sim.linSolveMethod=im; printf( "COMMAND LINE: -method(%i)   \n", im            ); } };
    //funcs["-dt"      ]={1,[&](const char** ss){ float dt;  sscanf( ss[0], "%f", &dt );  sim.dt=dt;             printf( "COMMAND LINE: -dt( dt: %f ) \n", sim.dt        ); } };
    //funcs["-bmix"    ]={1,[&](const char** ss){ int istart; float bmix;  sscanf( ss[0], "%i,%f", &istart, &bmix ); W.sim.mixer.b_end=bmix; W.sim.mixer.istart=istart; printf( "COMMAND LINE: -bmix( istart:%i bmix: %f ) \n", W.sim.mixer.istart, W.sim.mixer.b_end );    } };
    //funcs["-fix"     ]={1,[&](const char** ss){ int n =  readlist( ss[0], W.fixPoints); printf("COMMAND LINE: -fix[%i]{%s}\n", n, ss[0] );  } };
    //funcs["-nsolve"  ]={1,[&](const char** ss){ int nsolv; sscanf( ss[0], "%i", &nsolv ); printf( "COMMAND LINE: -nsolve(%i) \n", nsolv ); W.sim_f.nSolverIters=nsolv; W.sim.nSolverIters=nsolv;  } };    
    funcs["-parabola"] = {0, [&](const char**){ 
        printf("funcs[-parabola]: Parabola Extrude Test: "); 
        truss.clear();
        // --- Create truss by extruding Parabola (e.g. to build magnetic nozzle for nuclear spacecraft)
        // Parabola2Mesh(truss,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  true ); // Wireframe
        // Parabola2Mesh(truss,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  false ); // Solid (crashes due to normal calculation in MeshBuilder.h)
        Parabola_ExtrudedWire( truss, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.5, 0.1 );
        printf("Parabola Extrude Test: "); truss.printSizes();
    }};
    funcs["-skelet"]   = {0, [&](const char**){ 
        printf("funcs[-skelet]: Truss from Skelet Test:\n");
        skelet.clear();
        std::vector<double> node_sizes;
        createSkelet(skelet, node_sizes);
        trussFromSkelet(skelet, node_sizes, bb, truss);
        printf("Truss from Skelet Test:\n");
        printf("  skelet.printSizes(): "); skelet.printSizes();
        printf("  truss.printSizes():  "); truss.printSizes();
        skelet.write_obj("skelet.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.write_obj("truss.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges | ObjMask::Tris));
    
    }};
    funcs["-blocks"]   = {0, [&](const char**){ 
        printf("funcs[-blocks]: Manual Construction Blocks Test:\n");
        testConstructionBlocks(bb, truss);
        printf("Manual Construction Blocks Test: "); truss.printSizes();
    }};

    //if( argc<=1 ){  funcs["-parabola"].func(0); }
    //if( argc<=1 ){  funcs["-blocks"].func(0); }
    if( argc<=1 ){  funcs["-skelet"].func(0); }


    process_args( argc, argv, funcs );

	//app = new ConstructionBlockApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
