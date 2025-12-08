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

namespace SpaceCrafting {

// Global pointer to the editor's simulator used by SpaceCraftControl
static SpaceCraftSimulator* gEditorSimulator = nullptr;

void SpaceCraftControl(double dt){
    if(gEditorSimulator){
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

    // === Editor state
    EDIT_MODE edit_mode = EDIT_MODE::component;
    int picked = -1;
    Vec3d mouse_hray;
    Vec3d mouse_ray0;
    int picked_block = -1;
    bool bSmoothLines = true;
    bool bWireframe = true;

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
        SpaceCraftDynamicsApp::draw(); // Call parent's draw first
        
        // Additional editor-specific drawing
        if(edit_mode == EDIT_MODE::component) {
            if(compGui) compGui->draw();
        } else if(edit_mode == EDIT_MODE::edge) {
            if(girderGui) girderGui->draw();
        } else if(edit_mode == EDIT_MODE::vertex) {
            if(plateGui) plateGui->draw();
        }
        
        // Draw picking interface
        if(picked >= 0) {
            glColor3f(0,1.0,0);
            if(_sim) drawPicked(*theSpaceCraft, picked);
        }
    }

    virtual void drawHUD() override {
        SpaceCraftDynamicsApp::drawHUD();
        gui.draw();
        
        // Draw editor mode info
        char str[256];
        const char* mode_name = "unknown";
        switch(edit_mode) {
            case EDIT_MODE::component: mode_name = "Component"; break;
            case EDIT_MODE::vertex:    mode_name = "Vertex";    break;
            case EDIT_MODE::edge:      mode_name = "Edge";      break;
        }
        sprintf(str, "Edit Mode: %s", mode_name);
        Draw::drawText(str, fontTex, fontSizeDef, {WIDTH,HEIGHT-20});
    }

    virtual void eventHandling(const SDL_Event& event) override {
        SpaceCraftDynamicsApp::eventHandling(event);
        gui.onEvent(mouseX, mouseY, event);
        
        switch(event.type) {
        case SDL_MOUSEBUTTONDOWN:
            if(event.button.button == SDL_BUTTON_LEFT) {
                mouse_ray0 = (Vec3d)cam.pos;
                mouse_hray = (Vec3d)cam.rot.c;  // Forward direction from camera                
                if(edit_mode == EDIT_MODE::component && _mesh) {
                    // Configure picker ray before picking and ensure picker backend is valid
                    picker.ray0 = mouse_ray0;
                    picker.hray = mouse_hray;
                    if(picker.picker) {
                        picker.pick(); // Use picker's built-in pick function
                    }
                }
            }
            break;
        case SDL_KEYDOWN:
            switch(event.key.keysym.sym) {
                case SDLK_m: picker.switch_mode(); break;
                case SDLK_l: onSelectLuaShipScript.GUIcallback(lstLuaFiles); break;
            }
            break;
        }
    }

    virtual void keyStateHandling(const Uint8 *keys) override {
        SpaceCraftDynamicsApp::keyStateHandling(keys);
        
        if(keys[SDL_SCANCODE_1]) edit_mode = EDIT_MODE::component;
        if(keys[SDL_SCANCODE_2]) edit_mode = EDIT_MODE::vertex;
        if(keys[SDL_SCANCODE_3]) edit_mode = EDIT_MODE::edge;
        
        if(keys[SDL_SCANCODE_W]) bWireframe = !bWireframe;
        if(keys[SDL_SCANCODE_L]) bSmoothLines = !bSmoothLines;
    }

    void reloadShip(const char* fname) {
        if(simulator) {
            simulator->reloadShip(fname);
            // Initialize simulators using the mesh built from the Lua ship
            simulator->initSimulators();

            // Setup sliders and dynamic control, analogously to spaceCraftEditor.cpp and spaceCraftDynamics.cpp
            for( Ring* o : theSpaceCraft->rings ){ o->updateSlidersPaths( true, true, simulator->sim.points ); }
            sliders2edgeverts( *theSpaceCraft, simulator->sim );
            simulator->sim.user_update = SpaceCraftControl;
        }
    }
};

} // namespace SpaceCrafting

// ===================== MAIN



int main(int argc, char *argv[]) {
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
