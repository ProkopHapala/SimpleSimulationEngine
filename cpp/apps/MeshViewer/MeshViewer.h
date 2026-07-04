
#ifndef MeshViewer_h
#define MeshViewer_h

/// @file MeshViewer.h
/// @brief SDL2/OpenGL mesh visualizer with folder browser, animation playback, and SVG export.
///
/// Built on the engine's SDL2OGL infrastructure (AppSDL2OGL_3D, Draw2D, Draw3D, GUI).
/// Uses fixed-function OpenGL for rendering. Mouse/camera handling inherited from AppSDL2OGL_3D.
/// Text rendering via Draw2D::drawText (HUD) and Draw3D::drawInt (3D billboarded labels).
///
/// Features:
/// - **Folder browser**: TreeView with single-click navigation. Mesh files (.obj, .obj+, .npz) load on click.
/// - **Rendering**: solid (lit, with polygon offset), wireframe (derived from triangles if no edges), points.
/// - **Index labels**: billboarded vertex/edge/face numbers via Draw3D::drawInt (constant screen size).
/// - **Face normals**: red arrows from face centroids.
/// - **GUI**: CheckBoxList for display toggles, frame slider for multi-snapshot files.
/// - **Camera**: right-drag orbit (base class), wheel zoom, arrow keys rotate, numpad preset views.
/// - **SVG export**: faces/edges/points/normals/labels with backface culling, opacity, projection axis.

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <algorithm>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Camera.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"

#include "MeshFileFormats.h"  // now in common/geometry/ (global include path)

class MeshViewerApp : public AppSDL2OGL_3D { public:

    GUI gui;
    TreeView* treeView = nullptr;
    std::string currentDir;
    bool browseMode = true;

    MeshAnimation anim;
    std::string loadedFile;
    bool autoPlay = false;
    int playTimer = 0;
    int playInterval = 30;

    bool showSolid = true;
    bool showWire = true;
    bool showPoints = true;
    bool showVertNums = false;
    bool showEdgeNums = false;
    bool showFaceNums = false;
    bool showNormals = false;

    // GUI controls
    CheckBoxList* cboxes = nullptr;
    GUIPanel* frameSlider = nullptr;

    // SVG export options
    SVGExportOpts svgOpts;

    Vec3f camTarget = {0, 0, 0};
    float camFitDist = 20.0f;

    MeshViewerApp(int& id, int WIDTH_, int HEIGHT_) : AppSDL2OGL_3D(id, WIDTH_, HEIGHT_) {
        fontTex       = makeTexture    ("common_resources/dejvu_sans_mono_RGBA_inv.bmp");
        GUI_fontTex   = makeTextureHard("common_resources/dejvu_sans_mono_RGBA_pix.bmp");
        Draw::fontTex = fontTex;
        camDist = 20.0f;
        cam.zoom = 10.0f;
        cam.aspect = (float)WIDTH / HEIGHT;
        cam.persp = false; // orthographic
        cam.zmin = -10000.0f;
        cam.zmax =  10000.0f;

        // Checkbox panel for display options
        cboxes = new CheckBoxList(320, HEIGHT - 20*fontSizeDef, 580, fontSizeDef*2);
        cboxes->addBox("solid",  &showSolid);
        cboxes->addBox("wire",   &showWire);
        cboxes->addBox("points", &showPoints);
        cboxes->addBox("vnums",  &showVertNums);
        cboxes->addBox("enums",  &showEdgeNums);
        cboxes->addBox("fnums",  &showFaceNums);
        cboxes->addBox("normals", &showNormals);
        gui.addPanel(cboxes);
    }

    void setDir(const std::string& dir) {
        currentDir = dir;
        if (!treeView) {
            treeView = new TreeView("DirView", 5, 5, 300, HEIGHT / (fontSizeDef * 2) - 2);
            gui.addPanel(treeView);
            treeView->onActivate = [this](TreeViewTree* node) {
                if (!node->content.bDir && isMeshFile(node->content.caption))
                    loadFile(node->content.path);
            };
        }
        // clear old tree
        for (auto* br : treeView->root.branches) delete br;
        treeView->root.branches.clear();
        dir2tree(treeView->root, dir.c_str(), "");
        treeView->updateLines();
        treeView->iSelected = 0;
        treeView->redraw = true;
    }

    void loadFile(const std::string& path) {
        for (auto* s : anim.snapshots) delete (MeshSnapshot*)s;
        anim.snapshots.clear();
        anim.current = 0;
        if (!loadMeshFile(path, anim)) { printf("Failed to load %s\n", path.c_str()); return; }
        loadedFile = path;
        fitCameraToMesh();
        CMesh* snap = anim.currentSnap();
        if (snap) printf("Loaded %s: %d snapshots | verts=%d edges=%d faces=%d\n", path.c_str(), anim.size(), snap->nvert, snap->nedge, snap->ntri);

        // Setup frame slider if multiple snapshots
        if (frameSlider) { gui.panels.erase(std::remove(gui.panels.begin(), gui.panels.end(), frameSlider), gui.panels.end()); delete frameSlider; frameSlider = nullptr; }
        if (anim.size() > 1) {
            frameSlider = new GUIPanel("frame", 320, HEIGHT - 22*fontSizeDef, 580, HEIGHT - 20*fontSizeDef, true, true, true, true, false);
            frameSlider->setRange(0, anim.size() - 1);
            frameSlider->setValue(0);
            frameSlider->command = [this](GUIAbstractPanel* p){
                GUIPanel* gp = (GUIPanel*)p;
                anim.current = gp->getIntVal();
                fitCameraToMesh();
            };
            gui.addPanel(frameSlider);
        }
    }

    void fitCameraToMesh() {
        CMesh* snap = anim.currentSnap();
        if (!snap || snap->nvert == 0) return;
        Vec3d c = cmesh_center(snap);
        double r = cmesh_boundingRadius(snap, c);
        if (r < 1e-6) r = 1.0;
        camTarget = {(float)c.x, (float)c.y, (float)c.z};
        camFitDist = (float)(r * 3.0 + 1.0);
        cam.zoom = (float)(r * 1.2); // ortho: half-height of view = 1.2 * radius
        cam.pos = camTarget + Vec3f{0, 0, -camFitDist};
        cam.rot.setOne();
        qCamera.setOne();
        cam.zmin = -camFitDist * 2;
        cam.zmax =  camFitDist * 2;
    }

    virtual void draw() override {
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        if (anim.size() > 0) drawMesh();
    }

    virtual void drawHUD() override {
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        gui.draw();
        if (anim.size() > 0) drawOverlay();
        else {
            Draw::setRGB(0x404040);
            const char* help = "Up/Down: navigate  Enter: open dir / load file  Backspace: parent  Esc: quit";
            Draw2D::drawText(help, strlen(help), {5, fontSizeDef * 2}, 0.0f, GUI_fontTex, fontSizeDef);
        }
    }

    void drawMesh() {
        if (autoPlay) {
            playTimer++;
            if (playTimer >= playInterval) { playTimer = 0; anim.next(); }
        }

        camera();
        glEnable(GL_DEPTH_TEST);

        CMesh* snap = anim.currentSnap();
        if (!snap) return;

        if (showSolid && snap->ntri > 0) {
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0f, 1.0f);
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
            GLfloat lightPos[] = {0.5f, 0.5f, 1.0f, 0.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
            GLfloat lightDiff[] = {0.7f, 0.7f, 0.65f, 1.0f};
            glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff);
            GLfloat lightAmb[] = {0.4f, 0.4f, 0.45f, 1.0f};
            glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmb);

            glBegin(GL_TRIANGLES);
            for (int i = 0; i < snap->ntri; i++) {
                Vec3i iv = snap->tris[i];
                Vec3d a = snap->verts[iv.a], b = snap->verts[iv.b], c = snap->verts[iv.c];
                Vec3d nor; nor.set_cross(b - a, c - a); nor.normalize();
                Draw3D::normal(nor);
                Draw3D::vertex(a); Draw3D::vertex(b); Draw3D::vertex(c);
            }
            glEnd();
            glDisable(GL_LIGHTING);
            glDisable(GL_POLYGON_OFFSET_FILL);
        }

        if (showWire) {
            glDisable(GL_DEPTH_TEST); // always draw edges on top
            glColor3f(0.1f, 0.4f, 0.1f);
            glBegin(GL_LINES);
            if (snap->nedge > 0) {
                for (int i = 0; i < snap->nedge; i++) {
                    Draw3D::vertex(snap->verts[snap->edges[i].x]);
                    Draw3D::vertex(snap->verts[snap->edges[i].y]);
                }
            } else if (snap->ntri > 0) {
                for (int i = 0; i < snap->ntri; i++) {
                    Vec3i iv = snap->tris[i];
                    Draw3D::vertex(snap->verts[iv.a]); Draw3D::vertex(snap->verts[iv.b]);
                    Draw3D::vertex(snap->verts[iv.b]); Draw3D::vertex(snap->verts[iv.c]);
                    Draw3D::vertex(snap->verts[iv.c]); Draw3D::vertex(snap->verts[iv.a]);
                }
            }
            glEnd();
            glEnable(GL_DEPTH_TEST);
        }

        if (showPoints) {
            glDisable(GL_DEPTH_TEST);
            glPointSize(4.0f);
            glColor3f(0.8f, 0.4f, 0.0f);
            glBegin(GL_POINTS);
            for (int i = 0; i < snap->nvert; i++) Draw3D::vertex(snap->verts[i]);
            glEnd();
            glEnable(GL_DEPTH_TEST);
        }

        if (showNormals && snap->ntri > 0) {
            glDisable(GL_DEPTH_TEST);
            glColor3f(0.8f, 0.0f, 0.0f);
            glBegin(GL_LINES);
            for (int i = 0; i < snap->ntri; i++) {
                Vec3i iv = snap->tris[i];
                Vec3d a = snap->verts[iv.a], b = snap->verts[iv.b], c = snap->verts[iv.c];
                Vec3d cen = (a + b + c) * (1.0/3.0);
                Vec3d nor; nor.set_cross(b - a, c - a); nor.normalize();
                double ns = cam.zoom * 0.1; // normals scale with zoom (world-space, not text)
                Draw3D::vertex(cen);
                Draw3D::vertex(cen + nor * ns);
            }
            glEnd();
            glEnable(GL_DEPTH_TEST);
        }

        // Index labels (drawn in 3D, billboarded)
        glDisable(GL_DEPTH_TEST);
        if (showVertNums) {
            Draw::setRGB(0x0000A0);
            Draw3D::pointLabels(snap->nvert, snap->verts, fontTex, 0.02);
        }
        if (showEdgeNums) {
            Draw::setRGB(0xA00000);
            if (snap->nedge > 0) {
                for (int i = 0; i < snap->nedge; i++) {
                    Vec3d mid = (snap->verts[snap->edges[i].x] + snap->verts[snap->edges[i].y]) * 0.5;
                    Draw3D::drawInt(mid, i, fontTex, 0.02);
                }
            } else if (snap->ntri > 0) {
                int ei = 0;
                for (int i = 0; i < snap->ntri; i++) {
                    Vec3i iv = snap->tris[i];
                    Vec3d mid;
                    mid = (snap->verts[iv.a] + snap->verts[iv.b]) * 0.5; Draw3D::drawInt(mid, ei++, fontTex, 0.02);
                    mid = (snap->verts[iv.b] + snap->verts[iv.c]) * 0.5; Draw3D::drawInt(mid, ei++, fontTex, 0.02);
                    mid = (snap->verts[iv.c] + snap->verts[iv.a]) * 0.5; Draw3D::drawInt(mid, ei++, fontTex, 0.02);
                }
            }
        }
        if (showFaceNums) {
            Draw::setRGB(0x00A000);
            for (int i = 0; i < snap->ntri; i++) {
                Vec3i iv = snap->tris[i];
                Vec3d centroid = (snap->verts[iv.a] + snap->verts[iv.b] + snap->verts[iv.c]) * (1.0/3.0);
                Draw3D::drawInt(centroid, i, fontTex, 0.02);
            }
        }
        glEnable(GL_DEPTH_TEST);
    }

    void drawOverlay() {
        Draw::setRGB(0x202060);
        std::string info = loadedFile;
        if (anim.size() > 1) info += "  [Frame " + std::to_string(anim.current + 1) + "/" + std::to_string(anim.size()) + "]";
        Draw2D::drawText(info.c_str(), info.length(), {320, HEIGHT - fontSizeDef * 2}, 0.0f, GUI_fontTex, fontSizeDef);

        CMesh* snap = anim.currentSnap();
        if (snap) {
            Draw::setRGB(0x404080);
            char buf[256];
            snprintf(buf, sizeof(buf), "verts: %d  edges: %d  faces: %d", snap->nvert, snap->nedge, snap->ntri);
            Draw2D::drawText(buf, strlen(buf), {320, HEIGHT - fontSizeDef * 4}, 0.0f, GUI_fontTex, fontSizeDef);
        }

        Draw::setRGB(0x404040);
        const char* help = "RMB:orbit Arrows:rotate Wheel:zoom Space:play [:prev ]:next W:wire S:solid P:pts F:fit R:reload X:exportSVG";
        Draw2D::drawText(help, strlen(help), {320, fontSizeDef * 2}, 0.0f, GUI_fontTex, fontSizeDef);

        Draw::setRGB(0x505080);
        const char* viewHelp = "Numpad: 0:front 1:back 2:top 3:bottom 4:left 5:right 6:default";
        Draw2D::drawText(viewHelp, strlen(viewHelp), {320, fontSizeDef * 10}, 0.0f, GUI_fontTex, fontSizeDef);

        Draw::setRGB(0x505080);
        char idxBuf[256];
        snprintf(idxBuf, sizeof(idxBuf), "Labels: V=%s E=%s F=%s  (7:vert 8:edge 9:face)",
            showVertNums ? "on" : "off", showEdgeNums ? "on" : "off", showFaceNums ? "on" : "off");
        Draw2D::drawText(idxBuf, strlen(idxBuf), {320, fontSizeDef * 4}, 0.0f, GUI_fontTex, fontSizeDef);

        // SVG options line
        Draw::setRGB(0x505080);
        char svgBuf[256];
        snprintf(svgBuf, sizeof(svgBuf), "SVG: faces=%s op=%.2f bkf=%s ax=%c%s%s vnum=%s enum=%s nrm=%s",
            svgOpts.showFaces ? "on" : "off",
            svgOpts.faceOpacity,
            svgOpts.showBackfaces ? "on" : "off",
            svgOpts.projAxis == 0 ? 'X' : svgOpts.projAxis == 1 ? 'Y' : 'Z',
            svgOpts.showEdges ? " E" : "",
            svgOpts.showPoints ? " P" : "",
            svgOpts.showVertNums ? "on" : "off",
            svgOpts.showEdgeNums ? "on" : "off",
            svgOpts.showNormals ? "on" : "off");
        Draw2D::drawText(svgBuf, strlen(svgBuf), {320, fontSizeDef * 6}, 0.0f, GUI_fontTex, fontSizeDef);

        const char* svgHelp = "X:export 1:faces 2:opacity 3:backface 4:edges 5:axis 6:pts N:vnums M:enums B:normals";
        Draw::setRGB(0x606060);
        Draw2D::drawText(svgHelp, strlen(svgHelp), {320, fontSizeDef * 8}, 0.0f, GUI_fontTex, fontSizeDef);

        if (autoPlay) {
            Draw::setRGB(0x40E040);
            const char* play = "PLAYING";
            Draw2D::drawText(play, strlen(play), {WIDTH - 80, HEIGHT - fontSizeDef * 2}, 0.0f, GUI_fontTex, fontSizeDef);
        }
    }

    virtual void camera() override {
        qCamera.toMatrix(cam.rot);
        cam.pos = camTarget + cam.rot.c * (-camFitDist);
        Cam::ortho(cam, false);
    }

    virtual void eventHandling(const SDL_Event& event) override {
        gui.onEvent(mouseX, mouseY, event);
        if (event.type == SDL_KEYDOWN) {
            handleBrowseKey(event.key.keysym.sym);
            handleViewKey(event.key.keysym.sym);
        }
        if (event.type == SDL_MOUSEWHEEL) {
            if (event.wheel.y > 0) cam.zoom *= 0.9f;
            else if (event.wheel.y < 0) cam.zoom *= 1.1f;
        }
        AppSDL2OGL_3D::eventHandling(event);
    }

    void handleBrowseKey(SDL_Keycode key) {
        if (!treeView) return;
        treeView->updateLines();
        int n = treeView->lines.size();
        if (key == SDLK_UP)        { if (treeView->iSelected > 0) treeView->iSelected--; treeView->redraw = true; }
        else if (key == SDLK_DOWN) { if (treeView->iSelected < n - 1) treeView->iSelected++; treeView->redraw = true; }
        else if (key == SDLK_RETURN || key == SDLK_KP_ENTER) {
            if (n > 0 && treeView->iSelected < n) {
                TreeViewTree* sel = treeView->lines[treeView->iSelected];
                if (sel->content.bDir) { sel->content.open ^= true; treeView->redraw = true; }
                else if (isMeshFile(sel->content.caption)) loadFile(sel->content.path);
            }
        }
        else if (key == SDLK_BACKSPACE) {
            size_t slash = currentDir.find_last_of('/');
            if (slash > 0) setDir(currentDir.substr(0, slash));
        }
    }

    void handleViewKey(SDL_Keycode key) {
        if (key == SDLK_SPACE)  { autoPlay = !autoPlay; printf("Auto-play: %s\n", autoPlay ? "ON" : "OFF"); return; }
        // Arrow keys handled by base class keyStateHandling for rotation
        if (key == SDLK_w) { showWire = !showWire; return; }
        if (key == SDLK_s) { showSolid = !showSolid; return; }
        if (key == SDLK_p) { showPoints = !showPoints; return; }
        if (key == SDLK_7) { showVertNums = !showVertNums; cboxes->syncRead(); cboxes->redraw=true; return; }
        if (key == SDLK_8) { showEdgeNums = !showEdgeNums; cboxes->syncRead(); cboxes->redraw=true; return; }
        if (key == SDLK_9) { showFaceNums = !showFaceNums; cboxes->syncRead(); cboxes->redraw=true; return; }
        if (key == SDLK_LEFTBRACKET)  { anim.prev(); fitCameraToMesh(); if(frameSlider){frameSlider->setValue(anim.current);} printf("Frame %d/%d\n", anim.current+1, anim.size()); return; }
        if (key == SDLK_RIGHTBRACKET) { anim.next(); fitCameraToMesh(); if(frameSlider){frameSlider->setValue(anim.current);} printf("Frame %d/%d\n", anim.current+1, anim.size()); return; }
        if (key == SDLK_f) { fitCameraToMesh(); return; }
        if (key == SDLK_r) { if (!loadedFile.empty()) loadFile(loadedFile); return; }

        // Preset camera views
        if (key == SDLK_KP_0 || key == SDLK_0) { qCamera = Quat4fIdentity; printf("View: front\n"); return; }
        if (key == SDLK_KP_1)                  { qCamera = Quat4fFront;    printf("View: back\n");  return; }
        if (key == SDLK_KP_2)                  { qCamera = Quat4fTop;      printf("View: top\n");   return; }
        if (key == SDLK_KP_3)                  { qCamera = Quat4fBotton;   printf("View: bottom\n");return; }
        if (key == SDLK_KP_4)                  { qCamera = Quat4fLeft;     printf("View: left\n");  return; }
        if (key == SDLK_KP_5)                  { qCamera = Quat4fRight;    printf("View: right\n"); return; }
        if (key == SDLK_KP_6)                  { qCamera = Quat4fBack;     printf("View: default\n");return; }

        // SVG export controls
        if (key == SDLK_x) { exportSVG(); return; }
        if (key == SDLK_1) { svgOpts.showFaces = !svgOpts.showFaces; printf("SVG faces: %s\n", svgOpts.showFaces ? "on" : "off"); return; }
        if (key == SDLK_2) {
            svgOpts.faceOpacity += 0.25; if (svgOpts.faceOpacity > 1.0) svgOpts.faceOpacity = 0.0;
            printf("SVG opacity: %.2f\n", svgOpts.faceOpacity); return;
        }
        if (key == SDLK_3) { svgOpts.showBackfaces = !svgOpts.showBackfaces; printf("SVG backfaces: %s\n", svgOpts.showBackfaces ? "on" : "off"); return; }
        if (key == SDLK_4) { svgOpts.showEdges = !svgOpts.showEdges; printf("SVG edges: %s\n", svgOpts.showEdges ? "on" : "off"); return; }
        if (key == SDLK_5) { svgOpts.projAxis = (svgOpts.projAxis + 1) % 3; printf("SVG proj axis: %d\n", svgOpts.projAxis); return; }
        if (key == SDLK_6) { svgOpts.showPoints = !svgOpts.showPoints; printf("SVG points: %s\n", svgOpts.showPoints ? "on" : "off"); return; }
        if (key == SDLK_n) { svgOpts.showVertNums = !svgOpts.showVertNums; printf("SVG vert nums: %s\n", svgOpts.showVertNums ? "on" : "off"); return; }
        if (key == SDLK_m) { svgOpts.showEdgeNums = !svgOpts.showEdgeNums; printf("SVG edge nums: %s\n", svgOpts.showEdgeNums ? "on" : "off"); return; }
        if (key == SDLK_b) { svgOpts.showNormals = !svgOpts.showNormals; printf("SVG normals: %s\n", svgOpts.showNormals ? "on" : "off"); return; }
    }

    void exportSVG() {
        CMesh* snap = anim.currentSnap();
        if (!snap || snap->nvert == 0) { printf("No mesh loaded\n"); return; }
        // Derive output path from loaded file name
        std::string outPath = "mesh_export.svg";
        if (!loadedFile.empty()) {
            size_t dot = loadedFile.find_last_of('.');
            size_t slash = loadedFile.find_last_of('/');
            std::string base = (slash != std::string::npos) ? loadedFile.substr(slash + 1) : loadedFile;
            dot = base.find_last_of('.');
            if (dot != std::string::npos) base = base.substr(0, dot);
            if (anim.size() > 1) base += "_frame" + std::to_string(anim.current);
            outPath = base + ".svg";
        }
        writeSVG(outPath, *snap, svgOpts);
    }

    virtual void keyStateHandling(const Uint8* keys) override {
        float moveSpeed = 0.5f;
        if (keys[SDL_SCANCODE_A]) camTarget.add_mul(cam.rot.a, -moveSpeed);
        if (keys[SDL_SCANCODE_D]) camTarget.add_mul(cam.rot.a,  moveSpeed);
        if (keys[SDL_SCANCODE_Q]) camFitDist *= 1.02f;
        if (keys[SDL_SCANCODE_E]) camFitDist *= 0.98f;
    }
};

#endif
