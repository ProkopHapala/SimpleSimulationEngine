#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "globals.h"

#include "SpaceCraftDynamicsApp.h"
#include "EditSpaceCraft.h"
#include "SpaceCraftGUI.h"
#include "spaceCraftEditorUtils.h"
#include "IO_utils.h"
#include "Tree.h"
#include "LuaHelpers.h"
#include "argparse.h"
#include "Solids.h"
#include "BucketUtils.h"
#include "SpaceCraftDraw.h"

namespace SpaceCrafting {

// Global pointer to the editor's simulator used by SpaceCraftControl
static SpaceCraftSimulator* gEditorSimulator = nullptr;

void SpaceCraftControl(double dt){
    if(gEditorSimulator){
        // if(verbosity>0){
        //     int nSliders = (int)theSpaceCraft->sliders.size();
        //     double pathCur0 = (nSliders>0) ? theSpaceCraft->sliders[0]->path.cur : 0.0;
        //     double evc0 = (gEditorSimulator->sim.nEdgeVertBonds>0) ? gEditorSimulator->sim.edgeVertBonds[0].c : 0.0;
        //     printf("SpaceCraftEditorNew::SpaceCraftControl(dt=%g): wheel_speed=(%g,%g,%g) nSliders=%d path0.cur=%g ev0.c=%g\n",  dt,gEditorSimulator->wheel_speed.x, gEditorSimulator->wheel_speed.y, gEditorSimulator->wheel_speed.z, nSliders, pathCur0, evc0);
        // }
        applySliders2sim( *theSpaceCraft, gEditorSimulator->sim, (double*)&gEditorSimulator->wheel_speed );
    }
}

class SpaceCraftEditorNew : public SpaceCraftDynamicsApp {
public:
    // === Editor-specific GUI elements
    BoundGUI*  compGui   = nullptr;
    GirderGUI* girderGui = nullptr;
    PlateGUI*  plateGui  = nullptr;
    PickerUI   picker;
    DropDownList* lstLuaFiles = nullptr;
    CheckBoxList* dbgBoxes = nullptr;

    // === Editor state
    EDIT_MODE edit_mode = EDIT_MODE::component;
    int picked = -1;
    Vec3d mouse_hray;
    Vec3d mouse_ray0;
    int picked_block = -1;
    bool bSmoothLines = true;
    bool bWireframe = true;
    bool bVertexNumbers = false;
    bool bDebugSliders  = false;
    bool bHover         = false;
    bool bBBoxDebug     = false;
    bool bPolyWireframe = false;
    bool bShowHelp      = false;
    int  hoverPicked    = -1;
    bool bDragging      = false;
    int  dragVertex     = -1;
    Vec3d p_debug{ 100.0, 0.0, 0.0 };
    char str_tmp[8096];

    class OnSelectLuaShipScript : public GUIEventCallback { 
    public:
        SpaceCraftEditorNew* editor;
        OnSelectLuaShipScript(SpaceCraftEditorNew* editor_) : editor(editor_) {}
        virtual int GUIcallback(GUIAbstractPanel* caller) override {
            char str[256];
            ((DropDownList*)caller)->selectedToStr(str + sprintf(str,"data/"));
            editor->reloadShip(str);
            return 0;
        }
    };
    OnSelectLuaShipScript onSelectLuaShipScript;

    // === Constructor
    SpaceCraftEditorNew(int& id, int WIDTH_, int HEIGHT_, int argc, char *argv[]) : 
        SpaceCraftDynamicsApp(id, WIDTH_, HEIGHT_),
        onSelectLuaShipScript(this)
    {
        // Initialize GUI elements
        fontTex       = makeTexture    ("common_resources/dejvu_sans_mono_RGBA_inv.bmp");
        GUI_fontTex   = makeTextureHard("common_resources/dejvu_sans_mono_RGBA_pix.bmp");
        Draw::fontTex = fontTex;

        plateGui  = (PlateGUI* )gui.addPanel( new PlateGUI ( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );
        girderGui = (GirderGUI*)gui.addPanel( new GirderGUI( WIDTH-105, 5, WIDTH-5, fontSizeDef*2+2) );
        
        lstLuaFiles = new DropDownList("lua files", 20, HEIGHT-100, 200, 5); 
        gui.addPanel(lstLuaFiles);
        lstLuaFiles->setCommand([this](GUIAbstractPanel* panel){ onSelectLuaShipScript.GUIcallback(panel); });

        // --- Debug toggle checkboxes
        dbgBoxes = new CheckBoxList( 5, 5, 160, fontSizeDef*2 );
        dbgBoxes->addBox( "Vertex nums [N]", &bVertexNumbers );
        dbgBoxes->addBox( "Debug sliders [B]", &bDebugSliders  );
        dbgBoxes->addBox( "Hover [H]",        &bHover         );
        dbgBoxes->addBox( "BBox debug [V]",   &bBBoxDebug     );
        dbgBoxes->addBox( "Wireframe [P]",    &bPolyWireframe );
        gui.addPanel( dbgBoxes );

        // Load ship files list using Tree functionality
        TreeView* tvDir = new TreeView("DirView", 20, HEIGHT-400, 200, 20);
        gui.addPanel(tvDir);
        dir2tree(tvDir->root, "data");
        tvDir->updateLines();

        // Initialize display lists
        ogl_asteroide = glGenLists(1);
        glNewList(ogl_asteroide, GL_COMPILE);
        drawAsteroide(32, 50, 0.1, false);
        glEndList();

        ogl_geoSphere = glGenLists(1);
        glNewList(ogl_geoSphere, GL_COMPILE);
        drawAsteroide(16, 0, 0.0, true);
        glEndList();

        // Load lua files list
        listDirContaining("data", ".lua", lstLuaFiles->labels);

        // Setup view parameters
        VIEW_DEPTH = 10000.0;
        zoom = 1000.0;

        // Initialize simulator and spacecraft
        simulator = new SpaceCraftSimulator();
        // Bind simulator to base app and connect picker after _sim is valid
        bindSimulators(simulator);
        if(_sim) picker.picker = _sim;
        picker.Rmax = 10.0;

        // Expose simulator to SpaceCraftControl callback
        gEditorSimulator = simulator;

        // Initialize workshop (materials)
        init_workshop();

        // Initialize Lua and optionally load default ship when no CLI args
        initSpaceCraftingLua();
        if(argc <= 1) {
            reloadShip("data/ship_ICF_marksman_2.lua");
        }
    }

    virtual void draw() override {
        glLineWidth(1.0); // reset line width (may be stale from previous frame)
        if(bPolyWireframe){ glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); } else { glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); }
        SpaceCraftDynamicsApp::draw(); // Call parent's draw first

        // --- Update spatial buckets for picking (if initialized)
        if(_sim){
            if(_sim->pointBBs.ncell>0) updatePointBBs( _sim->pointBBs, _sim->BBs, _sim->points,            true );
            if(_sim->edgeBBs .ncell>0) updateEdgeBBs ( _sim->edgeBBs,  _sim->BBs, _sim->bonds, _sim->points, false );
        }

        // --- Debug sliders
        if(bDebugSliders && _sim) debug_sliders();

        // --- BBox debug visualization
        if(bBBoxDebug && _sim){
            picked_block = _sim->pick_BBox( picker.ray0, picker.hray, 10000.0, 1 );
            renderPickedBBox( picked_block, *_sim );
        }

        // --- Vertex numbering overlay
        if(bVertexNumbers && _sim){
            glDisable(GL_DEPTH_TEST);
            glColor3f(0.0,0.0,0.0);
            for(int i=0; i<_sim->nPoint; i++){ Draw3D::drawInt( _sim->points[i].f, i, fontTex, 0.02 ); }
            glEnable(GL_DEPTH_TEST);
        }

        // --- Draw sliders
        if(_sim){
            glLineWidth(5.0);
            glColor3f(0.0,0.5,1.0); drawSliderBonds( *_sim );
            glColor3f(1.0,0.0,1.0); drawSliders    ( *theSpaceCraft, *_sim );
            glLineWidth(3.0);
            glColor3f(0.0,0.5,0.0); drawSliderPaths( *theSpaceCraft, *_sim );
            glLineWidth(1.0);
        }

        // --- Picking visualization
        if(_sim){
            picker.hray = (Vec3d)(cam.rot.c);
            picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
            glLineWidth(5.0);
            if     (picker.edit_mode == EDIT_MODE::vertex){ if( picker.picked>=0 ){ Vec3d p = *(Vec3d*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawPointCross( p, 10.0 );                              } }
            else if(picker.edit_mode == EDIT_MODE::edge  ){ if( picker.picked>=0 ){ Vec2i b = *(Vec2i*)picker.getPickedObject(); glColor3f(0.0,1.0,0.0); Draw3D::drawLine      ( _sim->points[b.x].f, _sim->points[b.y].f ); } }
            glLineWidth(1.0);
        }

        // --- Hover highlight
        if(bHover && hoverPicked>=0 && _sim){
            glColor3f(1.0,0.5,0.0);
            if     (picker.edit_mode == EDIT_MODE::vertex){ Draw3D::drawPointCross( _sim->points[hoverPicked].f, 8.0 ); }
            else if(picker.edit_mode == EDIT_MODE::edge  ){ int2 b = _sim->bonds[hoverPicked]; Draw3D::drawLine( _sim->points[b.x].f, _sim->points[b.y].f ); }
        }

        // --- Drag feedback
        if(bDragging && dragVertex >= 0 && _sim){
            glColor3f(1.0,0.0,0.0);
            Draw3D::drawPointCross( _sim->points[dragVertex].f, 12.0 );
        }

        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        if( (picker.picked>=0) && (edit_mode==EDIT_MODE::component) ){
            glColor3f(0,1.0,0);
            drawPicked( *theSpaceCraft, picked );
        }
    }

    virtual void drawHUD() override {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);  // reset wireframe so text/GUI renders correctly
        SpaceCraftDynamicsApp::drawHUD();

        if(dbgBoxes) dbgBoxes->syncRead();  // sync checkboxes with keyboard-toggled bools

        // --- Status line
        if(_sim){
            sprintf(str_tmp, "time=%10.5f[s] mass=%g cog(%g,%g,%g) vcog(%g,%g,%g) L(%g,%g,%g) torq(%g,%g,%g) |F|=%g \n", _sim->time, _sim->mass, _sim->cog.x,_sim->cog.y,_sim->cog.z, _sim->vcog.x,_sim->vcog.y,_sim->vcog.z, _sim->L.x,_sim->L.y,_sim->L.z, _sim->torq.x,_sim->torq.y,_sim->torq.z, _sim->F_residual );
            Draw::drawText( str_tmp, fontTex, fontSizeDef,  {WIDTH,HEIGHT-20}  );
        }

        // --- Mode + debug state line
        const char* modeNames[] = {"?","vertex","edge","?","component"};
        int modeIdx = (int)picker.edit_mode;
        const char* modeName = (modeIdx>=0 && modeIdx<=4) ? modeNames[modeIdx] : "?";
        sprintf(str_tmp, "mode=%s picked=%i hover=%s nums=%s sliders=%s bbox=%s wire=%s \n",
            modeName, picker.picked,
            bHover?"on":"off", bVertexNumbers?"on":"off", bDebugSliders?"on":"off",
            bBBoxDebug?"on":"off", bPolyWireframe?"on":"off" );
        Draw::drawText( str_tmp, fontTex, fontSizeDef,  {WIDTH,HEIGHT-40}  );

        // --- Help overlay
        if(bShowHelp){
            sprintf(str_tmp,
                "=== Key Shortcuts (F1 to close) ===\n"
                "M     - cycle pick mode (vertex/edge/component)\n"
                "N     - toggle vertex numbering\n"
                "B     - toggle slider debug viz\n"
                "H     - toggle hover highlight\n"
                "V     - toggle bbox debug viz\n"
                "P     - toggle wireframe/fill\n"
                "L     - reload Lua ship script\n"
                "Space - pause/resume simulation\n"
                "KP_0  - clear velocities\n"
                "[ ]   - cycle bbox blocks\n"
                "LMB   - pick / drag vertex (in vertex mode)\n"
                "RMB   - rotate camera\n"
                "Wheel - zoom\n"
                "F1    - toggle this help\n"
            );
            Draw::drawText( str_tmp, fontTex, fontSizeDef, {200, HEIGHT-40} );
        }

        gui.draw();
    }

    virtual void eventHandling(const SDL_Event& event) override {
        SpaceCraftDynamicsApp::eventHandling(event);
        // Note: gui.onEvent already called by base class

        switch(event.type) {
        case SDL_KEYDOWN:
            switch(event.key.keysym.sym) {
                case SDLK_LEFTBRACKET  : if(_sim) circ_inc(picked_block, _sim->edgeBBs.ncell ); break;
                case SDLK_RIGHTBRACKET : if(_sim) circ_dec(picked_block, _sim->edgeBBs.ncell ); break;
                case SDLK_m: picker.switch_mode(); break;
                case SDLK_F1: bShowHelp = !bShowHelp; break;
                case SDLK_n: bVertexNumbers = !bVertexNumbers; break;
                case SDLK_b: bDebugSliders  = !bDebugSliders;  break;
                case SDLK_h: bHover         = !bHover;         break;
                case SDLK_v: bBBoxDebug     = !bBBoxDebug;     break;
                case SDLK_p: bPolyWireframe  = !bPolyWireframe;  if(bPolyWireframe){ glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); } else { glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); } break;
                case SDLK_l: onSelectLuaShipScript.GUIcallback(lstLuaFiles); break;
                case SDLK_KP_0: if(_sim) _sim->cleanVel(); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            if(event.button.button == SDL_BUTTON_LEFT) {
                picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
                picker.hray = (Vec3d)(cam.rot.c);
                if(picker.picker) {
                    picker.pick();
                    if(picker.edit_mode == EDIT_MODE::vertex && picker.picked >= 0){
                        bDragging  = true;
                        dragVertex = picker.picked;
                    }
                }
            }
            break;
        case SDL_MOUSEBUTTONUP:
            if(event.button.button == SDL_BUTTON_LEFT) {
                bDragging  = false;
                dragVertex = -1;
            }
            break;
        case SDL_MOUSEMOTION:
            if(bDragging && dragVertex >= 0 && _sim){
                Vec3d hray = (Vec3d)(cam.rot.c);
                Vec3d ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
                Vec3d p0   = _sim->points[dragVertex].f;
                double t   = (p0 - ray0).dot(hray);
                _sim->points[dragVertex].f = ray0 + hray * t;
            }
            if(bHover && !bDragging && _sim){
                picker.hray = (Vec3d)(cam.rot.c);
                picker.ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
                if     (picker.edit_mode == EDIT_MODE::vertex){ hoverPicked = (_sim->pointBBs.ncell>0 && _sim->BBs) ? pick_point_bucket( _sim->pointBBs, _sim->BBs, _sim->points, picker.ray0, picker.hray, picker.Rmax ) : _sim->pick_point_brute(picker.ray0, picker.hray, picker.Rmax); }
                else if(picker.edit_mode == EDIT_MODE::edge  ){ hoverPicked = (_sim->edgeBBs .ncell>0 && _sim->BBs) ? pick_bond_bucket ( _sim->edgeBBs,  _sim->BBs, _sim->bonds, _sim->points, picker.ray0, picker.hray, picker.Rmax ) : _sim->pick_bond_brute (picker.ray0, picker.hray, picker.Rmax); }
                else { hoverPicked = -1; }
            } else { hoverPicked = -1; }
            break;
        }
    }

    virtual void keyStateHandling(const Uint8 *keys) override {
        SpaceCraftDynamicsApp::keyStateHandling(keys);

        if(keys[SDL_SCANCODE_1]) edit_mode = EDIT_MODE::component;
        if(keys[SDL_SCANCODE_2]) edit_mode = EDIT_MODE::vertex;
        if(keys[SDL_SCANCODE_3]) edit_mode = EDIT_MODE::edge;
    }

    void reloadShip(const char* fname) {
        if(simulator) {
            simulator->reloadShip(fname);
            // Initialize simulators using the mesh built from the Lua ship
            simulator->initSimulators();

            // Setup sliders and dynamic control
            for( Ring* o : theSpaceCraft->rings ){ o->updateSlidersPaths( true, true, simulator->sim.points ); }
            sliders2edgeverts( *theSpaceCraft, simulator->sim );
            simulator->sim.user_update = SpaceCraftControl;
            simulator->sim.updateInveriants(false);
        }
    }

    void debug_sliders(){
        printf( "============= debug_sliders() \n" );
        glLineWidth(5.0);
        Draw3D::drawPointCross( p_debug, 30.0 );
        uint32_t colors[4]{ 0xFF0000FF, 0xFF00FF00, 0xFFFF0000, 0xFFFFFF00};
        for( Slider* o: theSpaceCraft->sliders){
            o->print(); printf(" - boundTo:");
            o->boundTo->print();
            int side = o->boundTo->nearSide(p_debug);
            printf(" - side: %i \n", side );
            glBegin(GL_LINES);
            Draw::setRGB( colors[side] );
            for(int i=0; i<5; i++){
                Vec3d p;
                o->boundTo->pointAlong( 0.2*i+0.1, side, &p );
                glVertex3d( p.x, p.y, p.z );
                glVertex3d( p_debug.x, p_debug.y, p_debug.z  );
            }
            glEnd();
            for(int side=0; side<4; side++ ){
                glBegin(GL_LINE_STRIP);
                Draw::setRGB( colors[side] );
                for(int i=0; i<5; i++){
                    Vec3d p;
                    o->boundTo->pointAlong( 0.2*i+0.1, side, &p );
                    glVertex3d( p.x, p.y, p.z );
                }
                glEnd();
            }
        }
        glLineWidth(1.0);
    }

    void renderPickedBBox( int picked_block, TrussDynamics_d& sim ){
        if(picked_block>=0){
            glPointSize(5);
            glLineWidth(3);
            Quat8d bb = sim.BBs[picked_block];
            glColor3f( 0.0,1.0,0.0 ); Draw3D::drawBBox( bb.lo.f, bb.hi.f );
            glColor3f( 0.0,1.0,1.0 ); renderEdgeBox ( picked_block, sim.edgeBBs,  sim.bonds, sim.points );
            glColor3f( 1.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointChunks,         sim.points );
            glColor3f( 0.0,0.0,1.0 ); renderPointBox( picked_block, sim.pointBBs,            sim.points );
        }
        glLineWidth(1);
    }
};

} // namespace SpaceCrafting

// ===================== MAIN



int main(int argc, char *argv[]) {
    verbosity = 2;

    SDL_DisplayMode dm = initSDLOGL( 8 );
	int junk;
	SpaceCrafting::SpaceCraftEditorNew* app = new SpaceCrafting::SpaceCraftEditorNew( junk, dm.w-150, dm.h-100, argc, argv );

    // Handle command-line arguments (e.g., -s data/ship.lua)
    LambdaDict funcs;
    funcs["-s"] = {1, [&](const char** ss){ app->reloadShip( ss[0] ); }};
    process_args( argc, argv, funcs );

    app->loop(1000000);
    
    return 0;
}
