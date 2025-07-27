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
Mesh::ConstructionBlockToMeshBuilder cbm(&truss);

Mesh::Renderer trussView(truss);

const int GIRDER_TYPE = 1;
const int ROPE_TYPE   = 2;


SDF_Cylinder* sdf_cylinder; 
Vec3d pivot_point{ 5.0, 0.0, 3.0 };

Selection* sel1; 
Selection* sel2;


using namespace Mesh;

void trussFromSkelet(const Mesh::Builder2& skelet, double* node_sizes, BlockBuilder& bb, Mesh::Builder2& truss) {
    // skelet (Mesh::Builder2) -> BlockBuilder
    bb.clear();
    for (int i = 0; i < skelet.verts.size(); ++i) { bb.addBlock(skelet.verts[i].pos, node_sizes[i]);}
    for (const auto& edge : skelet.edges)         { if (edge.w == GIRDER_TYPE) { bb.connectBlocks(edge.x, edge.y); }}
    Mesh::ConstructionBlockToMeshBuilder cbm(&truss);
    // BlockBuilder -> truss (Mesh::Builder2)
    cbm.drawBlockBuilder(bb, 4);
    // ropes
    for (const auto& edge : skelet.edges) {
        // int rope ( Vec3d p0,  Vec3d p1, int n, int ropeType, int anchorType, double Rcolapse=0.1, double r=-1.0 );
        if (edge.w == ROPE_TYPE) { truss.rope(skelet.verts[edge.x].pos, skelet.verts[edge.y].pos, 4, -1,-1,   0.1,1.0 ); } // -1 for auto-segmentation
    }
}

void bridgeLineSelections( Vec3d a0, Vec3d a1, Vec3d b0, Vec3d b1, double r, bool bTris=false){ 
    // sdf_cylinder = new SDF_Cylinder(  a0, a1, r, false );
    // Selection* s1 = truss.curSelection; truss.selectVertsBySDF( *sdf_cylinder                     ); truss.nextSelection(); s1->print();
    // Selection* s2 = truss.curSelection; truss.selectVertsBySDF( SDF_Cylinder(  b0, b1, r, false ) ); s2->print(); 
    // sel1=s1; sel2=s2;    
    Selection* s1 = truss.curSelection; truss.selectVertsBySDF( SDF_Cylinder(  a0, a1, r, false ) ); truss.nextSelection( ); 
    Selection* s2 = truss.curSelection; truss.selectVertsBySDF( SDF_Cylinder(  b0, b1, r, false ) );
    sel1=s1; sel2=s2;
    int n1 = s1->vec.size();
    int n2 = s2->vec.size();
    int n = (n1<n2)?n1:n2;
    SDF_point2 costf{ (a0+b0)*0.5 };
    //printf("==== bridgeLineSelections() n1 %i n2 %i n %i\n", n1, n2, n);
    truss.sortVertsBy( s1->vec, costf );    // s1->print();
    truss.sortVertsBy( s2->vec, costf );    // s2->print();
    truss.bridgeTriPatch( n, s1->vec.data(), s2->vec.data(), r, Quat4i{-1,-1,-1,-1}, bTris );
    
} // select edges by capsula

void brideFormSkeleton( Vec2i e1, Vec2i e2, Vec3d* ps, double* szs, double r, double h, const Vec3d* up_=0, bool bTris=false){ 
    Vec3d p0 = ps[e1.x]; Vec3d p1 = ps[e1.y];
    Vec3d p2 = ps[e2.x]; Vec3d p3 = ps[e2.y];
    Vec3d d1 = p1 - p0; d1.normalize();
    Vec3d d2 = p3 - p2; d2.normalize();
    Vec3d up; if(up_==0){ up.set_cross(d1,d2); up.normalize(); }else{ up = *up_; }
    Vec2d sz1{ szs[e1.x], szs[e1.y] };
    Vec2d sz2{ szs[e2.x], szs[e2.y] };
    // printf("brideFormSkeleton() e1(%i,%i) e2(%i,%i)\n", e1.x, e1.y, e2.x, e2.y);
    // printf("brideFormSkeleton() p0(%g,%g,%g) p1(%g,%g,%g) p2(%g,%g,%g) p3(%g,%g,%g)\n", p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z, p3.x, p3.y, p3.z);
    // printf("brideFormSkeleton() d1(%g,%g,%g) d2(%g,%g,%g) up(%g,%g,%g) r %g h %g\n", d1.x, d1.y, d1.z, d2.x, d2.y, d2.z, up.x, up.y, up.z);
    bridgeLineSelections( 
        p0+d2*(sz1.x)+up*(sz1.x*h),  p1+d2*(sz1.y)+up*(sz1.y*h), 
        p2+d1*(sz2.x)+up*(sz2.x*h),  p3+d1*(sz2.y)+up*(sz2.y*h), r, bTris );
} // select edges by capsula


// ====================== Class Definitions

class ConstructionBlockApp : public AppSDL2OGL_3D { public:

    int perFrame = 10;

	bool bSmoothLines = 1;
	bool bWireframe   = 1;

    // View control properties
    bool bViewSkelet       = false;
    bool bViewBlockBuilder = false;
    bool bViewTruss        = true;
    bool bViewPivotPoint   = true;
    bool bViewSelection    = true;
    bool bViewAxis         = true;

    
    // bool bViewEdges        = true;
    // bool bViewTris         = true;
    // bool bViewFaces        = true;
    // bool bViewFaceNormals  = false;
    // bool bViewPointLabels  = false;
    // bool bViewFaceLabels   = false;
    // bool bViewTriLabels    = false;


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

    // --- Single block testing
    // block.Ls=Vec3d{1.1,1.0,0.9};
    // for(int i=0;i<6;i++){
    //     block.faces[i].typ=1;
    // }
    // truss.printSizes();

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


    trussView.bViewPointLabels  = true;
    //trussView.bViewEdgeLabels   = true;
    //trussView.bViewVertNormals  = true;
    trussView.bViewFaceLabels   = true;
    trussView.bViewFaceNormals  = true;
    //trussView.bViewTriLabels    = true;
    bViewAxis = false;
    
    viewControls->addBox("Block Builder", &bViewBlockBuilder);
    viewControls->addBox("Pivot Point",   &bViewPivotPoint);
    viewControls->addBox("Selection",     &bViewSelection);
    viewControls->addBox("Skelet",        &bViewSkelet);
    viewControls->addBox("Truss",         &bViewTruss);
    viewControls->addBox("Edges",         &trussView.bViewEdges);
    viewControls->addBox("Faces",         &trussView.bViewFaces);
    viewControls->addBox("Tris",          &trussView.bViewTris);
    viewControls->addBox("Face Normals",  &trussView.bViewFaceNormals);
    viewControls->addBox("Vert Normals",  &trussView.bViewVertNormals);

    viewControls->addBox("Point Labels",  &trussView.bViewPointLabels);
    viewControls->addBox("Edge Labels",   &trussView.bViewEdgeLabels);
    viewControls->addBox("Face Labels",   &trussView.bViewFaceLabels);
    viewControls->addBox("Tri Labels",    &trussView.bViewTriLabels);

    viewControls->addBox("Axis",          &bViewAxis);

    gui.addPanel(viewControls);

    //mp= new MultiPanel( "Edit", gx.x0, ylay.x0, gx.x1, 0,-13); 

    // context menu for right mouse button
    contextMenu = new MultiPanel( "Context Menu", 0, 0, fontSizeDef*20, fontSizeDef*2, -1 );
    contextMenu->hideOnCommand = true;
    //contextMenu->addPanel( "select along line",  {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->ffl.print_nonbonded();   return 0; }; 

    // contextMenu->addButton( "selectVertsAlongPolyline", [&](GUIAbstractPanel* p){ 
    //     printf( "selectVertsAlongPolyline \n" );
    //     printf( "select along line BEFORE \n" ); truss.printSelectedVerts();
    //     truss.selectVertsAlongPolyline( 0.1, true ); 
    //     printf( "select along line AFTER \n" ); truss.printSelectedVerts();
    //     return 0; 
    // } );

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

    if(bViewTruss) {
        glEnable(GL_LIGHTING);
        trussView.draw();
        trussView.drawLabels();
        glLineWidth(1.0);
        glColor3f(0.0,0.7,0.0);
        if(bDragging && iDraggingButton==SDL_BUTTON_LEFT ){ drawMuseSelectionBox(); }
        glDisable(GL_DEPTH_TEST);
        if(bViewSelection) {
            //trussView.drawSelection(ipick);
            //printf("ConstructionBlockApp::draw() bViewSelection @sel1=%p @sel2=%p\n", sel1, sel2);
            if(sel1){ glColor3f(1.0,1.0,0.0); trussView.drawSelection(ipick, sel1); }
            if(sel2){ glColor3f(0.0,1.0,1.0); trussView.drawSelection(ipick, sel2); }
        }
    }
    if(bViewPivotPoint) { Draw3D::drawPointCross( pivot_point, 0.5 ); }
    if(bViewAxis) { glLineWidth(5.0); Draw3D::drawAxis(10.0); glLineWidth(1.0); }



    if(sdf_cylinder) {
        //Draw3D::drawCapsula( Vec3d p0, Vec3d p1, float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped );
        //Draw3D::drawCapsula( sdf_cylinder->p0, sdf_cylinder->p0 + sdf_cylinder->hdir*sdf_cylinder->l, sdf_cylinder->r, sdf_cylinder->r, 0.0, 1.0, 0.1, 8, true );
        Draw3D::drawCapsula( sdf_cylinder->p0, sdf_cylinder->p0 + sdf_cylinder->hdir*sdf_cylinder->l, sdf_cylinder->r, sdf_cylinder->r, 0.0, 1.0, 0.1, 8, false, GL_LINE_STRIP );
    }
        
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
    // unbuffered
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    printf( "argc %i \n", argc );
    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	ConstructionBlockApp * app = new ConstructionBlockApp( junk , dm.w-150, dm.h-100, argc, argv );
    //app->bindSimulators( &W ); 

    LambdaDict funcs;

    funcs["-parabola"] = {0, [&](const char**){ 
        printf("funcs[-parabola]: Parabola Extrude Test: "); 
        truss.clear();
        // --- Create truss by extruding Parabola (e.g. to build magnetic nozzle for nuclear spacecraft)
        // Parabola2Mesh(truss,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  true ); // Wireframe
        // Parabola2Mesh(truss,{6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0,  false ); // Solid (crashes due to normal calculation in MeshBuilder.h)
        //Parabola_Wire( truss, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.0 );
        //Parabola_Wire_new( truss, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*1.5}, 10.0, 10.0, 0.0, true,  true  );

        //Parabola_Wire_new( truss, {8,5}, Vec2f{0.1,0.0}, Vec2f{1.0,M_PI*1.5}, 10.0, 10.0, 0.0, false, false );
        //Parabola_Wire_new( truss, {8,5}, Vec2f{0.1,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, true, true );
        
        //Parabola_Wire_new( truss, {3,6}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::DEFAULT_WIRE );
        
        //Parabola_Wire_new( truss, {3,8}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::STAR | WireFlags::DIAGONAL1_EDGES | WireFlags::ALTERNATE_DIAG );
        //Parabola_Wire_new( truss, {3,6}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::TRIMESH );

        Parabola_Wire_new( truss, {3,6}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::STAR  );
        //Parabola_Wire_new( truss, {3,16}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::TRIMESH );

        //Parabola_Wire( truss, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.5 );
        //Parabola_ExtrudedWire( truss, {6,10}, Vec2f{0.0,0.0}, Vec2f{1.0,M_PI*2-0.1}, 10.0, 10.0, 0.5, 0.1 );

        truss.build_edgesOfVerts();
        
        truss.selectRectEdge( Vec3dMin, Vec3dMax );
        truss.printSelection();
        pivot_point = Vec3d{0.0,0.0,+3.0};
        truss.normalsTowardPoint( truss.curSelection->vec.size(), truss.curSelection->vec.data(), pivot_point );

        truss.bevel( truss.curSelection->vec.size(), truss.curSelection->vec.data(), 0.1, 0.5, 1 );


        // select edge strip
        truss.curSelection->clear();
        printf("!!!!!!!!!! select edge strip\n");
        truss.build_edgesOfVerts();
        truss.selection_mode = (int)Mesh::Builder2::SelectionMode::edge;
        truss.selectEdgeStrip2( 15, 0.9, {10,-1,-1} );
        truss.printSelection(true);

        printf("Parabola Extrude Test: "); truss.printSizes();
    }};

    funcs["-panel"] = {0, [&](const char**){ 
        printf("funcs[-panel]: Panel Extrude Test: "); 
        truss.clear();
        //truss.panel( {0.0,0.0,0.0}, {100.0,0.0,0.0}, {0.0,100.0,0.0}, {100.0,100.0,0.0}, {4,3}, 10.0, Quat4i{0,0,0,0} );
        //Tube(truss, {4,16}, {-1.0,0.0}, {1.0,1.5*M_PI}, {10.0,10.0}, 10.0, 2.0, Quat4i{0,0,0,0} );
        //QuadPanel(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {80.0,0.0,0.0}, {0.0,100.0,0.0}, {120.0,120.0,0.0}, 10.0, Quat4i{0,0,0,0} );
        //QuadPanel(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 10.0, Quat4i{0,0,0,0} );
        //QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {80.0,0.0,0.0}, {0.0,100.0,0.0}, {120.0,120.0,0.0}, 10.0, 0b111,  Quat4i{0,0,0,0} );
        //QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 10.0, 0b111, Quat4i{0,0,0,0} );
        //QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 10.0, 0b1111, Quat4i{0,0,0,0} );
        //QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, {0.333333,0.333333,7.0}, 0b101010111, Quat4i{0,0,0,0} );
        //SlabTube(truss, {2,16}, {0.0,0.0}, {0.2,1.5*M_PI}, {10.0,10.0}, 10.0, {0.333333,0.333333,2.0}, 0b101010111, Quat4i{0,0,0,0} );
        //SlabTube(truss, {4,17}, {0.0,0.0}, {1.0,1.0}, {10.0,10.0}, 12.0, {0.0,0.5,2.0}, 0b101010111, Quat4i{0,0,0,0} );
        //SlabTube(truss, {4,9}, {0.0,0.0}, {1.0,1.0}, {10.0,10.0}, 24.0, {0.0,0.5,2.0}, 0b101010111, Quat4i{0,0,0,0} );

        //QuadSheet(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {80.0,0.0,0.0}, {0.0,100.0,0.0}, {120.0,120.0,0.0}, 0b1111,  Quat4i{0,0,0,0} );

        //QuadSheet(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 0b1011, Quat4i{0,0,0,0}, 3,10 );
        //QuadSheet(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 0b1011, Quat4i{0,0,0,0}, 3,15 );

        //TubeSheet(truss, {4,10}, {0.0,0.0}, {1.0,1.0}, {10.0,10.0}, 10.0, 0b1011, 0.5 );
        //TubeSheet(truss, {4,10}, {0.0,0.0}, {1.0,1.0}, {10.0,10.0}, 10.0, 0b1011, 0.0 );

        //TorusSheet(truss, {10,20}, {0.0,0.0}, {1.0,1.0}, {5.0,20.0}, 0b1011, 0.5 );    
        //TorusSheet(truss, {3,6}, {0.0,0.0}, {1.0,1.0}, {5.0,20.0}, 0b0011, 0.0 );
        //TorusSheet(truss, {3,6}, {-0.5,0.0}, {0.5,1.0}, {5.0,20.0}, 0b0011, 0.0 );
        TorusSheet(truss, {4,6}, {-0.125,0.0}, {0.875,1.0}, {5.0,20.0}, 0b0011, 0.0 );
        //TorusSheet(truss, {4,6}, {0.0,0.0}, {1.0,1.0}, {5.0,20.0}, 0b0011, 0.0 );

    }};

    funcs["-bevel"] = {0, [&](const char**){ 
        const int npoint      = 4;
        //const int npoint      = 7;
        Vec3d nodes[npoint] = {
            {0.0,    0.0,   0.0}, // 0
            {100.0,  0.0,   0.0}, // 1
            {-30.0, +80.0,  0.0}, // 2
            {-30.0, -80.0,  0.0}, // 3

            // {-30.0, +160.0,  0.0}, // 4
            // {150.0, +80.0,  0.0}, // 5
            // {150.0, -80.0,  0.0} // 6
        };
        //const int nedge       = 6;
        const int nedge       = 3;
        Vec2i edges[nedge] = {
            {0,1}, {0,2}, {0,3}
        //, {2,4}, {1,5}, {1,6} 
        };
        
        truss.add_verts(npoint, nodes);
        truss.add_edges(nedge, edges);
        truss.build_edgesOfVerts();

        //const int nbev = 6;
        const int nbev = 3;
        int ies[nbev] = {
            0,1,2
            //,3,4,5
        };

        truss.select_verts_of_edge( nbev, ies ); // select vertices of edges to curSelection
        truss.printSelection();
        pivot_point = {0.0, 0.0, 100.0};
        truss.normalsTowardPoint( truss.curSelection->vec.size(), truss.curSelection->vec.data(), pivot_point );

        truss.bevel( nbev, ies, 10.0, 10.0, 5.0);
    }};

    funcs["-skelet"]   = {1, [&](const char** ss){ 
        int bUseSkelet = (ss[0][0]=='T')||(ss[0][0]=='t')||(ss[0][0]=='1'); 
        printf("funcs[-skelet]: Truss from Skelet(%i):\n", bUseSkelet);
        skelet.clear();
        //std::vector<double> node_sizes;
        //createSkelet(skelet, node_sizes);
        const int N_NODES = 7;
        Vec3d node_positions[N_NODES] = {
            {  0.0,   0.0,   0.0}, // 0
            {  0.0,   0.0, -15.0}, // 1
            {  0.0,   0.0,  20.0}, // 2
            {  0.0, -10.0,   0.0}, // 3
            {  0.0,  10.0,   0.0}, // 4
            {-10.0,   0.0,   0.0}, // 5
            { 10.0,   0.0,   0.0}  // 6
        };
        double node_sizes[N_NODES] = {1.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.25};
        // Define girder edges
        const int N_GIRDER_EDGES = 6;
        Vec2i girder_edges[N_GIRDER_EDGES] = { {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6} };
        // Define rope edges
        const int N_ROPE_EDGES = 12;
        Vec2i rope_edges[N_ROPE_EDGES] = {
            {1, 3}, {1, 4}, {1, 5}, {1, 6},
            {2, 3}, {2, 4}, {2, 5}, {2, 6},
            {3, 5}, {5, 4}, {4, 6}, {6, 3}
        };
        if(bUseSkelet){
            printf("Using low-res skelet to create truss\n");
            bb.clear();
            skelet.addVerts( N_NODES, node_positions );
            skelet.ropes   ( N_GIRDER_EDGES, 1, girder_edges, GIRDER_TYPE );
            skelet.ropes   ( N_ROPE_EDGES,   1, rope_edges,   ROPE_TYPE   );
            printf("  skelet.printSizes(): "); skelet.printSizes();
            skelet.write_obj("skelet.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
            trussFromSkelet(skelet, node_sizes, bb, truss);
            //cbm.drawBlockBuilder(bb, 4);
        }else{
            printf("Using arrays directly to create truss\n");
            bb.addBlockNetwork( N_NODES, node_positions, node_sizes, N_GIRDER_EDGES, girder_edges );
            cbm.drawBlockBuilder(bb, 4);
            //truss.ropes( N_ROPE_EDGES, 4, rope_edges, ROPE_TYPE );
            // int ropes( int nv, Vec3d* vs, int ne, int nseg, const Vec2i* ends, int ropeType, int anchorType, double Rcolapse=0.1, double r=-1.0 );
            truss.ropes( N_NODES, node_positions, N_ROPE_EDGES, 4, rope_edges, ROPE_TYPE, GIRDER_TYPE, 0.1, 1.0 );

        }
        printf("  truss.printSizes():  "); truss.printSizes();
        truss .write_obj("truss.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges | ObjMask::Tris));
    
    }};

    funcs["-blocks"]   = {0, [&](const char**){ 
        printf("funcs[-blocks]: Manual Construction Blocks Test:\n");
        //testConstructionBlocks(bb, truss);
        const int N_NODES = 5;
        Vec3d node_positions[N_NODES] = {
            {0.0,   0.0,   0.0}, 
            {0.0,   0.0, -15.0}, 
            {0.0,   0.0,  20.0}, 
            {0.0, -10.0,   0.0}, 
            {0.0,  10.0,   0.0}
        };
        double node_sizes[N_NODES] = {1.0, 0.5, 0.5, 0.25, 0.25};
        const int N_GIRDER_EDGES = 4;
        Vec2i girder_edges[N_GIRDER_EDGES] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
        bb.addBlockNetwork(N_NODES, node_positions, node_sizes, N_GIRDER_EDGES, girder_edges );
        cbm.drawBlockBuilder(bb, 4);
        printf("Manual Construction Blocks Test: "); truss.printSizes();
    }};

    funcs["-extrude_octahedron"] = {0, [&](const char**){
        printf("--- Running test: extrude_octahedron\n");
        truss.addCMesh(Solids::Octahedron, false); // bFaces=false, we only want the wireframe initially
        Vec2i chs = truss.addFaces( Solids::Octahedron_nplanes, Solids::Octahedron_planes, Solids::Octahedron_planeVs, true );
        printf("Extruding chunk %i by 2.0 units...\n", chs.a);
        truss.extrudeFace(chs.a, 2.0);
    }};

    funcs["-oct_nodes"] = {0, [&](const char**){
        printf("funcs[-oct_nodes]: Manual Construction Blocks Test:\n");
        //testConstructionBlocks(bb, truss);
        const int nnodes = 5;
        Vec3d node_positions[nnodes] = {
            {0.0,   0.0,   0.0}, 
            {0.0,   0.0, -15.0}, 
            {0.0,   0.0,  20.0}, 
            {0.0, -10.0,   0.0}, 
            {0.0,  10.0,   0.0}
        };
        double node_sizes[nnodes] = {1.0, 0.5, 0.5, 0.25, 0.25};
        const int nedges = 4;
        Vec2i edges[nedges] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};
        Vec2i chs[nnodes];
        bool bUseSpecialPlanes=true;
        if(bUseSpecialPlanes){
            truss.facingNodes( Solids::Octahedron, nnodes, node_positions, chs, node_sizes, Solids::Octahedron_nplanes, Solids::Octahedron_planes, Solids::Octahedron_planeVs );
        }else{
            CMesh oct=(CMesh){Solids::Octahedron_nverts,Solids::Octahedron_nedges,Solids::Octahedron_ntris,Solids::Octahedron_nplanes, Solids::Octahedron_verts, Solids::Octahedron_edges, Solids::Octahedron_tris, Solids::Octahedron_planes, Solids::Octahedron_planeVs};
            truss.facingNodes( oct, nnodes, node_positions, chs, node_sizes );
        }
        //for(int i=0;i<nnodes;i++){ chs[i].y++; } // make range end exclusive
        truss.bridgeFacingPolygons( nedges, edges, node_positions, 4, chs );
        printf("  truss.printSizes():  "); truss.printSizes();
        truss.write_obj("truss.obj", ObjMask::Verts | ObjMask::Edges | ObjMask::Polygons );
        //printf("Extruding chunk %i by 2.0 units...\n", chs.a);
        //truss.extrudeFace(chs.a, 2.0);
    }};

    funcs["-cube_nodes"] = {0, [&](const char**){
        printf("funcs[-oct_nodes]: Manual Construction Blocks Test:\n");
        const int nnodes = 5;
        Vec3d node_ps[nnodes] = {
            {0.0,   0.0,   0.0}, 
            {0.0,   0.0, -15.0}, 
            {0.0,   0.0,  20.0}, 
            {0.0, -10.0,   0.0}, 
            {0.0,  10.0,   0.0}
        };
        double node_szs[nnodes] = {1.0, 0.5, 0.5, 0.25, 0.25};
        const int nedges = 4;
        const int nseg   = 5;
        Vec2i edges[nedges] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}};

        // const int nnodes = 2;
        // Vec3d node_ps[nnodes] = {
        //     {0.0,   0.0,   0.0}, 
        //     {0.0,   0.0,   3.0}
        // };
        // double node_szs[nnodes] = {1.0, 0.5};
        // const int nedges = 1;
        // const int nseg=1;
        // Vec2i edges[nedges] = {{0, 1}};

        Vec2i chs[nnodes];
        truss.facingNodes( Solids::Cube, nnodes, node_ps, chs, node_szs );
        truss.printFaces();
        truss.bridgeFacingPolygons( nedges, edges, node_ps, nseg, chs );
        //truss.bridgeFacingPolygons( nedges, edges, node_positions, nseg, chs, {-1,-1,-1,-1}, {0,0,0,0} );
        //printf("  truss.printSizes():  "); truss.printSizes();

        //brideFormSkeleton( {0,1}, {0,3}, node_ps, node_szs, 0.5, 1.0 );
        
        Vec3d up = {1,0,0};
        bool bTris=true;
        brideFormSkeleton( edges[0], edges[2], node_ps, node_szs, 0.5, 1.0, &up, bTris );
        brideFormSkeleton( edges[0], edges[3], node_ps, node_szs, 0.5, 1.0, &up, bTris );
        brideFormSkeleton( edges[1], edges[2], node_ps, node_szs, 0.5, 1.0, &up, bTris );
        brideFormSkeleton( edges[1], edges[3], node_ps, node_szs, 0.5, 1.0, &up, bTris );

        // const int nvs=3;
        // int ivs1[nvs] = {43,47,51};
        // int ivs2[nvs] = {90,94,98};
        // truss.triPatch( nvs, ivs1, ivs2, 0.5, {0,0,0,0} );

        //truss.bridgeVertStrips( nvs, ivs1, ivs2 );




        truss.write_obj("truss.obj", ObjMask::Verts | ObjMask::Edges | ObjMask::Polygons );
    }};

    //if( argc<=1 ){  funcs["-parabola"].func(0); }
    //if( argc<=1 ){  funcs["-blocks"].func(0); }
    if( argc<=1 ){  funcs["-skelet"].func(0); }

    process_args( argc, argv, funcs );

	//app = new ConstructionBlockApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
