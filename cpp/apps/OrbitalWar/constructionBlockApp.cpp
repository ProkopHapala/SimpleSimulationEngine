#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <array>
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
#include "KinematicSolver.h"
#include "MeshFileFormats.h"
#include "argparse.h"

int read_binary_int( const char* s, const char* label=0 ){
    int n=0;
    for(int i=0;i<32;i++){
        if(s[i]=='\0') break;
        if(s[i]=='1') n |= 1<<i;
    }
    if(label) printf("read_binary_int() %s: %i\n", label, n);
    return n;
}

Vec2i scanf_Vec2i( const char* s, const char* label=0 ){
    Vec2i v; sscanf(s, "%i,%i", &v.x, &v.y); 
    if(label) printf("scanf_Vec2i() %s: %i %i\n", label, v.x, v.y);
    return v;
}

Vec2f scanf_Vec2f( const char* s, const char* label=0 ){
    Vec2f v; sscanf(s, "%f,%f", &v.x, &v.y); 
    if(label) printf("scanf_Vec2f() %s: %f %f\n", label, v.x, v.y);
    return v;
}

Vec3f scanf_Vec3f( const char* s, const char* label=0 ){
    Vec3f v; sscanf(s, "%f,%f,%f", &v.x, &v.y, &v.z); 
    if(label) printf("scanf_Vec3f() %s: %f %f %f\n", label, v.x, v.y, v.z);
    return v;
}

bool bHeadless = false;

float scanf_float( const char* s, const char* label=0 ){
    float v = atof(s);
    if(label) printf("scanf_float() %s: %f\n", label, v);
    return v;
}


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

ConstructionBlockApp::ConstructionBlockApp( int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[] ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_, makeWindowTitle("ConstructionBlockApp", argc, argv).c_str() ) {
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

    funcs["-QuadSlab"] = {7, [&](const char** ss){ 
        printf("funcs[-QuadSlab]: Panel Extrude Test: \n"); 
        // QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, {0.333333,0.333333,7.0}, 0b101010111, Quat4i{0,0,0,0} );
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec3f p00   = scanf_Vec3f(ss[1], "p00");   
        Vec3f p01   = scanf_Vec3f(ss[2], "p01");   
        Vec3f p10   = scanf_Vec3f(ss[3], "p10");   
        Vec3f p11   = scanf_Vec3f(ss[4], "p11");   
        Vec3f up    = scanf_Vec3f(ss[5], "up");
        int dirMask = read_binary_int(ss[6], "dirMask"); 
        //double offset; sscanf(ss[5], "%lf", &offset);
        truss.clear();
        QuadSlab(truss, nuv, {0.0,0.0}, {1.0,1.0}, p00, p01, p10, p11, up, dirMask );
        //QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, {0.333333,0.333333,7.0}, 0b101010111, Quat4i{0,0,0,0} );
    }};


    funcs["-SlabTube"] = {6, [&](const char** ss){ 
        printf("funcs[-SlabTube]: Panel Extrude Test: \n"); 
        //SlabTube(truss, {2,16}, {0.0,0.0}, {0.2,1.5*M_PI}, {10.0,10.0}, 10.0, {0.333333,0.333333,2.0}, 0b101010111, Quat4i{0,0,0,0} );
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec2f UVmin = scanf_Vec2f(ss[1], "UVmin"); 
        Vec2f UVmax = scanf_Vec2f(ss[2], "UVmax"); 
        Vec3f RLs   = scanf_Vec3f(ss[3], "RLs");   
        Vec3f up    = scanf_Vec3f(ss[4], "up");
        int dirMask = read_binary_int(ss[5], "dirMask"); 
        //double offset; sscanf(ss[5], "%lf", &offset);
        truss.clear();
        //SlabTube(truss, nuv, UVmin, UVmax, RLs.xy(), RLs.z, up, dirMask );
        SlabTube(truss, {2,16}, {0.0,0.0}, {0.2,1.5*M_PI}, {10.0,10.0}, 10.0, {0.333333,0.333333,2.0}, 0b101010111, Quat4i{0,0,0,0} );
    }};

    funcs["-QuadSheet"] = {7, [&](const char** ss){ 
        printf("funcs[-QuadSheet]: Panel Extrude Test: \n"); 
        //QuadSheet(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 0b1011, Quat4i{0,0,0,0}, 3,15 );
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec3f p00   = scanf_Vec3f(ss[1], "p00");   
        Vec3f p01   = scanf_Vec3f(ss[2], "p01");   
        Vec3f p10   = scanf_Vec3f(ss[3], "p10");   
        Vec3f p11   = scanf_Vec3f(ss[4], "p11");    
        int dirMask = read_binary_int(ss[5], "dirMask"); 
        Vec2i iminmax = scanf_Vec2i(ss[6], "iminmax");
        //double offset; sscanf(ss[5], "%lf", &offset);
        truss.clear();
        QuadSheet(truss, nuv, Vec2f{0.0,0.0}, Vec2f{1.0,1.0}, p00, p01, p10, p11, dirMask, Quat4i{0,0,0,0}, iminmax.x, iminmax.y);
    }};

    funcs["-TubeSheet"] = {5, [&](const char** ss){ 
        printf("funcs[-TubeSheet]: Panel Extrude Test: \n"); 
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec2f UVmin = scanf_Vec2f(ss[1], "UVmin"); 
        Vec2f UVmax = scanf_Vec2f(ss[2], "UVmax"); 
        Vec3f RLs   = scanf_Vec3f(ss[3], "RLs");   
        int dirMask = read_binary_int(ss[4]); 
        //double offset; sscanf(ss[5], "%lf", &offset);
        truss.clear();
        TubeSheet(truss, nuv, UVmin, UVmax, RLs.xy(), RLs.z, dirMask, 0.0);
    }};

    funcs["-TorusSheet"] = {5, [&](const char** ss){ 
        printf("funcs[-TorusSheet]: Panel Extrude Test: \n"); 
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec2f UVmin = scanf_Vec2f(ss[1], "UVmin"); 
        Vec2f UVmax = scanf_Vec2f(ss[2], "UVmax"); 
        Vec2f Rs    = scanf_Vec2f(ss[3], "Rs");   
        int dirMask = read_binary_int(ss[4]); 
        truss.clear();
        TorusSheet(truss, nuv, UVmin, UVmax, Rs, dirMask, 0.0);
    }};

    funcs["-ParabolaSheet"] = {6, [&](const char** ss){ 
        printf("funcs[-ParabolaSheet]: Parabolic single-layer truss: \n"); 
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec2f UVmin = scanf_Vec2f(ss[1], "UVmin"); 
        Vec2f UVmax = scanf_Vec2f(ss[2], "UVmax"); 
        float R     = scanf_float(ss[3], "R");     
        float L     = scanf_float(ss[4], "L");     
        int dirMask = read_binary_int(ss[5]); 
        truss.clear();
        ParabolaSheet(truss, nuv, UVmin, UVmax, R, L, dirMask, 0.5);
        printf("ParabolaSheet: "); truss.printSizes();
        truss.write_obj("parabola_sheet.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("parabola_sheet.svg");
        // Multi-view: top, front, side, 45deg
        Mat3d rots[4] = { Mat3dIdentity, Mat3d{0,0,1, 0,1,0, -1,0,0}, Mat3d{1,0,0, 0,0,-1, 0,1,0}, Mat3d{0.7071,0,0.7071, 0,1,0, -0.7071,0,0.7071} };
        truss.exportSVGmultiView("parabola_sheet_views.svg", 4, rots);
    }};

    funcs["-ParabolaSlab"] = {7, [&](const char** ss){ 
        printf("funcs[-ParabolaSlab]: Parabolic double-layer slab: \n"); 
        Vec2i nuv   = scanf_Vec2i(ss[0], "nuv");   
        Vec2f UVmin = scanf_Vec2f(ss[1], "UVmin"); 
        Vec2f UVmax = scanf_Vec2f(ss[2], "UVmax"); 
        float R     = scanf_float(ss[3], "R");     
        float L     = scanf_float(ss[4], "L");     
        Vec3f up    = scanf_Vec3f(ss[5], "up");    
        int dirMask = read_binary_int(ss[6]); 
        truss.clear();
        ParabolaSlab_wrap(truss, nuv, UVmin, UVmax, R, L, up, dirMask, -0.5);
        printf("ParabolaSlab_wrap: "); truss.printSizes();
        truss.write_obj("parabola_slab.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("parabola_slab.svg");
        Mat3d rots[4] = { Mat3dIdentity, Mat3d{0,0,1, 0,1,0, -1,0,0}, Mat3d{1,0,0, 0,0,-1, 0,1,0}, Mat3d{0.7071,0,0.7071, 0,1,0, -0.7071,0,0.7071} };
        truss.exportSVGmultiView("parabola_slab_views.svg", 4, rots);
    }};

    funcs["-ParametricParabola"] = {6, [&](const char** ss){ 
        printf("funcs[-ParametricParabola]: Annular parabolic patch: \n"); 
        int   nTop    = atoi(ss[0]); 
        int   nBottom = atoi(ss[1]); 
        int   nRows   = atoi(ss[2]); 
        float R1      = scanf_float(ss[3], "R1"); 
        float R2      = scanf_float(ss[4], "R2"); 
        float L       = scanf_float(ss[5], "L");   
        truss.clear();
        ParametricParabolaPatch(truss, nTop, nBottom, nRows, R1, R2, L);
        printf("ParametricParabolaPatch: "); truss.printSizes();
        truss.write_obj("parabola_patch.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("parabola_patch.svg");
        Mat3d rots[4] = { Mat3dIdentity, Mat3d{0,0,1, 0,1,0, -1,0,0}, Mat3d{1,0,0, 0,0,-1, 0,1,0}, Mat3d{0.7071,0,0.7071, 0,1,0, -0.7071,0,0.7071} };
        truss.exportSVGmultiView("parabola_patch_views.svg", 4, rots);
    }};

    funcs["-TubeSheetBond"] = {0, [&](const char**){
        printf("funcs[-TubeSheetBond]: Recoil damper (hex outer + triangular inner) with parabolic dish\n");
        truss.clear();
        int dirMask = 0b1111;
        float twist = 0.0f;
        Quat4i stickTypes{1,0,0,0};

        // --- Outer hex tube ---
        Vec2i n1{2, 6};
        Vec2f UVmin{0,0}, UVmax{1,1};
        float R1 = 1.732f, L1 = 1.0f;
        Vec2f Rs1{R1, R1};
        int iv0_1 = truss.verts.size();
        TubeSheet(truss, n1, UVmin, UVmax, Rs1, L1, dirMask, twist, stickTypes);
        int nVerts1 = n1.x * n1.y;
        std::vector<int> sel1(nVerts1);
        for(int i=0; i<nVerts1; i++) sel1[i] = iv0_1 + i;

        // --- Second tube (shifted UV window) ---
        Vec2i n2{3, 6};
        float R2 = 1.0f, L2 = 2.0f;
        Vec2f Rs2{R2, R2};
        float offX = -0.5f, offY = 0.5f;
        float du = offX / (n2.x - 1.0f);
        float dv = offY / (float)n2.y;
        Vec2f UVmin2{du, dv}, UVmax2{1+du, 1+dv};
        int iv0_2 = truss.verts.size();
        TubeSheet(truss, n2, UVmin2, UVmax2, Rs2, L2, dirMask, twist, stickTypes);
        int nVerts2 = n2.x * n2.y;
        std::vector<int> sel2(nVerts2);
        for(int i=0; i<nVerts2; i++) sel2[i] = iv0_2 + i;

        // --- Bond-length analysis and connections ---
        double Rcut = 2.0, dR = 0.01;
        std::vector<std::vector<Vec2i>> groups;
        auto uniqLs = truss.mapBondLengthsFromVertexSelections(sel1, sel2, Rcut, dR, groups);
        if(!uniqLs.empty() && !groups.empty()){
            int k = std::max(0, std::min(0, (int)uniqLs.size()-1));
            truss.addEdgesFromPairs(groups[k], stickTypes.x);
        }

        // --- Inner triangular tube ---
        Vec2i nTri{10, 3};
        float clearance = 0.1f;
        float Rtri = std::max(0.0f, R2 - clearance);
        Vec2f RsTri{Rtri, Rtri};
        float dvTri = 1.0f / (2.0f * nTri.y) * 0.5f;
        float axshift = -0.9f;
        Vec2f UVminTri{axshift, dvTri}, UVmaxTri{1+axshift, 1+dvTri};
        TubeSheet(truss, nTri, UVminTri, UVmaxTri, RsTri, 20.0f, dirMask, 0.0f, stickTypes);

        // --- Parabolic pusher-plate / plasma nozzle dish ---
        int nRows=4, nInner=7, nOuter=25;
        float R1dish = 1.732f, R2dish = 1.732f*6, Ldish = 3.0f*6;
        float zShiftDish = L2*0.5f;
        int iv0_dish = truss.verts.size();
        ParametricParabolaPatch(truss, nOuter, nInner, nRows, R1dish, R2dish, Ldish);
        for(int i=iv0_dish; i<(int)truss.verts.size(); i++) truss.verts[i].pos.z += zShiftDish;

        printf("TubeSheetBond: "); truss.printSizes();
        truss.write_obj("tubesheet_bond.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("tubesheet_bond.svg");
        Mat3d rots[4] = { Mat3dIdentity, Mat3d{0,0,1, 0,1,0, -1,0,0}, Mat3d{1,0,0, 0,0,-1, 0,1,0}, Mat3d{0.7071,0,0.7071, 0,1,0, -0.7071,0,0.7071} };
        truss.exportSVGmultiView("tubesheet_bond_views.svg", 4, rots);
    }};

    funcs["-RopesVShape"] = {0, [&](const char**){
        printf("funcs[-RopesVShape]: Two polylines sharing a root vertex\n");
        truss.clear();
        int nseg = 10;
        Vec3d A{0,0,0}, B{0,10,0}, C{10,10,0};
        std::vector<int> vAB;
        for(int i=0; i<=nseg; i++){ float t=(float)i/nseg; vAB.push_back(truss.vert(A*(1-t)+B*t)); }
        for(int i=0; i<nseg; i++) truss.edge(vAB[i], vAB[i+1]);
        std::vector<int> vAC; vAC.push_back(vAB[0]);
        for(int i=1; i<=nseg; i++){ float t=(float)i/nseg; vAC.push_back(truss.vert(A*(1-t)+C*t)); }
        for(int i=0; i<nseg; i++) truss.edge(vAC[i], vAC[i+1]);
        printf("RopesVShape: "); truss.printSizes();
        truss.write_obj("ropes_vshape.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("ropes_vshape.svg");
    }};

    funcs["-RopesParallel"] = {0, [&](const char**){
        printf("funcs[-RopesParallel]: Two parallel polylines\n");
        truss.clear();
        int nseg = 10;
        Vec3d A1{-5,0,0}, B1{-5,10,0}, A2{5,0,0}, B2{5,10,0};
        for(auto& [P0,P1] : std::array<std::pair<Vec3d,Vec3d>,2>{{{A1,B1},{A2,B2}}}){
            std::vector<int> vs;
            for(int i=0; i<=nseg; i++){ float t=(float)i/nseg; vs.push_back(truss.vert(P0*(1-t)+P1*t)); }
            for(int i=0; i<nseg; i++) truss.edge(vs[i], vs[i+1]);
        }
        printf("RopesParallel: "); truss.printSizes();
        truss.write_obj("ropes_parallel.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("ropes_parallel.svg");
    }};

    funcs["-BridgeQuads"] = {0, [&](const char**){
        printf("funcs[-BridgeQuads]: Bridge two facing quads\n");
        truss.clear();
        float s = 2.0f;
        int v0=truss.vert({-s,-s,0}), v1=truss.vert({s,-s,0}), v2=truss.vert({s,s,0}), v3=truss.vert({-s,s,0});
        int v4=truss.vert({-s,-s,10}), v5=truss.vert({s,-s,10}), v6=truss.vert({s,s,10}), v7=truss.vert({-s,s,10});
        truss.bridge_quads({v0,v1,v2,v3}, {v4,v5,v6,v7}, 4, {0,1,2,3}, {1,1,1,1}, true);
        printf("BridgeQuads: "); truss.printSizes();
        truss.write_obj("bridge_quads.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("bridge_quads.svg");
    }};

    funcs["-BridgeQuadsTwisted"] = {0, [&](const char**){
        printf("funcs[-BridgeQuadsTwisted]: Bridge two quads with 45deg twist\n");
        truss.clear();
        float s = 2.0f; float ang = M_PI/4; float c=cos(ang), sn=sin(ang);
        int v0=truss.vert({-s,-s,0}), v1=truss.vert({s,-s,0}), v2=truss.vert({s,s,0}), v3=truss.vert({-s,s,0});
        int v4=truss.vert({-s*c-(-s)*sn, -s*sn+(-s)*c, 10});
        int v5=truss.vert({s*c-(-s)*sn, s*sn+(-s)*c, 10});
        int v6=truss.vert({s*c-s*sn, s*sn+s*c, 10});
        int v7=truss.vert({-s*c-s*sn, -s*sn+s*c, 10});
        truss.bridge_quads({v0,v1,v2,v3}, {v4,v5,v6,v7}, 6, {0,1,2,3}, {1,1,0,0}, true);
        printf("BridgeQuadsTwisted: "); truss.printSizes();
        truss.write_obj("bridge_quads_twisted.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("bridge_quads_twisted.svg");
    }};

    funcs["-ParametricQuadPatch"] = {0, [&](const char**){
        printf("funcs[-ParametricQuadPatch]: Variable-row triangulated patch\n");
        truss.clear();
        float s = 5.0f;
        Vec3d p00{-s,0,-s}, p01{-s,0,s}, p10{s,0,-s}, p11{s,0,s};
        ParametricQuadPatch(truss, 8, 20, 8, p00, p01, p10, p11);
        printf("ParametricQuadPatch: "); truss.printSizes();
        truss.write_obj("quad_patch.obj", (uint8_t)(ObjMask::Verts | ObjMask::Edges));
        truss.exportSVG("quad_patch.svg");
    }};

    funcs["-KinematicTest"] = {0, [&](const char**){
        printf("funcs[-KinematicTest]: damper demo — dish slides along triangular girder spine\n");
        truss.clear();

        int dirMask = 0b1111;
        Quat4i stickTypes{1,0,0,0};

        // --- Inner triangular girder (ship spine) — body 0, fixed ---
        // 3-sided cross-section, long along Z
        Vec2i nTri{10, 3};             // 10 axial, 3 circumferential (triangle)
        float R_girder = 0.9f;
        float L_girder = 20.0f;
        // Rotate triangle by half-segment via UV shift (like JS version)
        float dvTri = 1.0f / (2.0f * 3) * 0.5f;
        Vec2f UVminTri{0.0f, dvTri}, UVmaxTri{1.0f, 1.0f + dvTri};
        TubeSheet(truss, nTri, UVminTri, UVmaxTri, Vec2f{R_girder, R_girder}, L_girder, dirMask, 0.0f, stickTypes);
        int nTubeVerts = (int)truss.verts.size();

        // --- Parabolic pusher-plate dish — body 1, slides along girder ---
        int nRows=4, nInner=7, nOuter=25;
        float R1dish = 1.732f;         // inner hole radius (matches outer hex ring)
        float R2dish = 1.732f * 6;     // outer dish radius
        float Ldish  = 3.0f * 6;       // axial height
        float zShiftDish = 1.0f;       // start near back of girder
        int iv0_dish = truss.verts.size();
        ParametricParabolaPatch(truss, nOuter, nInner, nRows, R1dish, R2dish, Ldish);
        for(int i=iv0_dish; i<(int)truss.verts.size(); i++) truss.verts[i].pos.z += zShiftDish;

        printf("KinematicTest base mesh: "); truss.printSizes();

        // --- Setup kinematic system ---
        KinematicSystem ks;
        ks.poses.resize(2);
        ks.poses[0] = KinematicPose{Vec3d{0,0,0}, Mat3dIdentity, true};   // girder fixed
        ks.poses[1] = KinematicPose{Vec3d{0,0,0}, Mat3dIdentity, false};  // dish free

        // 3 slider constraints: paths derived from actual girder vertices
        // TubeSheet vertex layout (UV_slab_verts): idx[iy*n.x + ix], ix=axial(0..n.x-1), iy=circumferential(0..n.y-1)
        // Each rail = one circumferential side (fixed iy), varying ix along axial direction
        // This is exactly how Girder::sideToPath() works in SpaceCraftComponents.h
        int iv0_girder = 0;  // girder is first in truss
        std::vector<Vec3d>* pathVerts = new std::vector<Vec3d>();
        pathVerts->reserve(nTubeVerts);
        for(int i = 0; i < nTubeVerts; i++) pathVerts->push_back(truss.verts[i].pos);

        for(int s = 0; s < 3; s++){
            int iy = s;  // circumferential index = rail side
            std::vector<int>* pathIndices = new std::vector<int>(nTri.x);
            for(int ix = 0; ix < nTri.x; ix++) (*pathIndices)[ix] = iy * nTri.x + ix;

            KinematicConstraint c;
            c.type = ConstraintType::Slider;
            c.bodyA = 1;  // dish
            c.bodyB = -1; // path is world-fixed (girder is fixed)
            // Dish attachment point: at same (x,y) as the first vertex of this rail, z=0 in dish-local
            Vec3d railStart = truss.verts[iv0_girder + (*pathIndices)[0]].pos;
            c.lposA = Vec3d{railStart.x, railStart.y, 0};
            c.pathIndices = *pathIndices;
            c.pathClosed = false;
            c.pathVerts = pathVerts;
            c.sliderParamIndex = s;
            ks.constraints.push_back(c);

            printf("  rail %d: iy=%d, %d verts, start=(%.3f,%.3f,%.3f)\n", s, iy, (int)pathIndices->size(), railStart.x, railStart.y, railStart.z);
        }

        // All 3 sliders start at same t (dish at z=zShiftDish)
        ks.sliderParams.resize(3);
        double t0 = zShiftDish / L_girder;
        for(int s = 0; s < 3; s++) ks.sliderParams[s] = t0;

        printf("KinematicSystem: %d bodies, %d constraints, %d sliders, %d unknowns, %d equations\n",
               ks.nBodies(), (int)ks.constraints.size(), ks.nSliders(), ks.nUnknowns(), ks.nEquations());

        // --- Sweep: dish slides along girder rails ---
        int nSteps = 36;
        MeshAnimation anim;
        int nV = (int)truss.verts.size();
        int nE = (int)truss.edges.size();
        int nSliderEdges = (int)ks.constraints.size();
        // Rail path visualization: nPathPts-1 segments per rail (open path)
        int nRailVerts = 0, nRailEdges = 0;
        for(auto& c : ks.constraints){
            nRailVerts += (int)c.pathIndices.size();
            nRailEdges += (int)c.pathIndices.size() - 1;
        }
        double tStart = 0.05, tEnd = 0.85;
        ks.sweepSolveRange(tStart, tEnd, nSteps, [&](int step, const std::vector<KinematicPose>& poses){
            printf("  step %d: dish pos=(%.3f,%.3f,%.3f)\n", step, poses[1].pos.x, poses[1].pos.y, poses[1].pos.z);

            MeshSnapshot* snap = new MeshSnapshot();
            snap->alloc(nV + nSliderEdges * 2 + nRailVerts, nE + nSliderEdges + nRailEdges, 0);
            for(int i=0; i<nTubeVerts; i++) snap->verts[i] = truss.verts[i].pos;
            for(int i=iv0_dish; i<nV; i++) snap->verts[i] = localToWorld(poses[1], truss.verts[i].pos);
            for(int i=0; i<nE; i++){ snap->edges[i].x = truss.edges[i].lo.a; snap->edges[i].y = truss.edges[i].lo.b; }
            int vi = nV, ei = nE;
            // Rail paths (static lines showing where sliders should be)
            for(int s = 0; s < (int)ks.constraints.size(); s++){
                const auto& c = ks.constraints[s];
                int vStart = vi;
                for(int j = 0; j < (int)c.pathIndices.size(); j++){
                    snap->verts[vi] = (*c.pathVerts)[c.pathIndices[j]];
                    if(j > 0) snap->edges[ei++] = {vi - 1, vi};
                    vi++;
                }
                if(c.pathClosed) snap->edges[ei++] = {vi - 1, vStart};  // close the loop
            }
            // Slider connection lines: dish attachment point -> current path point
            for(int s = 0; s < (int)ks.constraints.size(); s++){
                const auto& c = ks.constraints[s];
                Vec3d dishPoint = localToWorld(poses[1], c.lposA);
                Vec3d pathPoint = bsplineInterpolate(ks.sliderParams[s], c.pathIndices, c.pathClosed, *c.pathVerts);
                snap->verts[vi]   = dishPoint;
                snap->verts[vi+1] = pathPoint;
                snap->edges[ei]   = {vi, vi+1};
                vi += 2; ei++;
            }
            anim.snapshots.push_back((CMesh*)snap);
        }, 200, 1e-4, true);

        writeOBJPlus("kinematic_anim.obj+", anim);
        printf("KinematicTest: done, wrote kinematic_anim.obj+ (%d frames)\n", anim.size());
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

    funcs["-headless"] = {0, [&](const char**){ 
        bHeadless = true;
        printf("funcs[-headless]: Headless mode enabled, will exit after export\n");
    }};

    //if( argc<=1 ){  funcs["-parabola"].func(0); }
    //if( argc<=1 ){  funcs["-blocks"].func(0); }
    if( argc<=1 ){  funcs["-skelet"].func(0); }

    process_args( argc, argv, funcs );

    if( bHeadless ){
        printf("Headless mode: exiting after export.\n");
        return 0;
    }

	//app = new ConstructionBlockApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
