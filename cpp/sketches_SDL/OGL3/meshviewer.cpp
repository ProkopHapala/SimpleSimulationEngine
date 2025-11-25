#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "Mat3.h"
#include "SDL_utils.h"
#include "Vec3.h"

#include "Mesh.h"

#include "AppSDL2OGL3.h"
#include "MeshRenderOGL3.h"
#include "SceneOGL3.h"
#include "ScreenSDL2OGL3.h"
#include "TextRendererOGL3.h"
#include "argparse.h"

static void checkGLError(const char *tag) {
  GLenum err;
  while ((err = glGetError()) != GL_NO_ERROR) {
    printf("GL_ERROR [%s]: 0x%04X\n", tag, err);
  }
}

class MeshViewerScene : public SceneOGL3 {
public:
  MeshRenderOGL3 renderer;
  OMesh mesh;
  Vec3f modelPos;

  TextRendererOGL3 text;
  bool labelZTest = false;

  void loadOBJ(const char *path) {
    mesh.fromFileOBJ(path);
    mesh.polygonsToTriangles(false);
    mesh.tris2normals(true);
    mesh.findEdges();
    renderer.uploadMesh_d(mesh.points.size(), mesh.triangles.size(),(int *)&mesh.triangles[0], (double *)&mesh.points[0], (double *)&mesh.normals[0]);
    Vec2i *evts = mesh.exportEdgeVs();
    renderer.uploadLines_d(mesh.points.size(), mesh.edges.size(), (int *)evts, (double *)&mesh.points[0]);
    delete[] evts;

    // Build labels: one label per vertex with its index as text
    text.clear();
    int n = mesh.points.size();
    for (int i = 0; i < n; i++) {
      char buf[16];
      sprintf(buf, "%d", i);
      Vec3f p = (Vec3f)mesh.points[i];
      //p.z += 0.1f; // small offset towards camera
      text.addLabel(p, buf);
    }
  }

  MeshViewerScene() {
    modelPos.set(0.0f, 0.0f, -30.0f);
    renderer.initDefaultShaders();
    renderer.setModelPos(modelPos);
    Mat3f ident;
    ident.setOne();
    renderer.setModelMat(ident);
    GLuint fontTex = makeTexture((char *)"common_resources/dejvu_sans_mono_RGBA_inv.bmp");
    text.init(fontTex, 0.5f, 0.5f);
    text.glyphOffset = 33;
    checkGLError("TextRendererOGL3::init");
  }

  virtual void draw(Camera &cam) override {
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    renderer.setModelPos(modelPos);
    renderer.draw(cam);

    if(!labelZTest) glDisable(GL_DEPTH_TEST);
    text.setModelPos(modelPos);
    text.draw(cam);
    if(!labelZTest) glEnable(GL_DEPTH_TEST);
    
    checkGLError("TextRendererOGL3::draw");
  }
};

class MeshViewerApp : public AppSDL2OGL3 {
public:
  MeshViewerScene *scene = nullptr;

  MeshViewerApp(int W, int H) : AppSDL2OGL3(W, H) {
    scene = new MeshViewerScene();
    screen->scenes.push_back(scene);
    // Basic camera setup
    screen->cam.pos.set(0.0f, 0.0f, 0.0f);
    screen->cam.rot.setOne();
    screen->cam.zmin = -1.0f;
    screen->cam.zmax = -1000.0f;
    screen->cam.zoom = 20.0f;
    screen->cam.aspect = W / (float)H;
    screen->cam.persp = true;

    // Orbit camera around the model center instead of rotating in-place
    screen->camLookAt = &scene->modelPos;
    screen->camDist = -30.0f; // roughly -scene->modelPos.z
  }

  virtual void eventHandling(const SDL_Event &event) override {
    // Keep default key handling (ESC, SPACE, QUIT)
    AppSDL2OGL3::eventHandling(event);

    // Add mouse-wheel zoom control
    if (event.type == SDL_MOUSEWHEEL && screen) {
      float zoomStep = 0.1f * screen->cam.zoom;
      if (event.wheel.y > 0) {
        screen->cam.zoom /= (1.0f + zoomStep);
      } else if (event.wheel.y < 0) {
        screen->cam.zoom *= (1.0f + zoomStep);
      }
    }
  }
};

int main(int argc, char *argv[]) {
  std::string objPath = "common_resources/turret.obj";
  int   winW      = -100;
  int   winH      = -100;
  float labelSize = 0.07f;
  float initZoom  = 20.0f;
  bool  doAutoFit = false;
  float lineWidth = 1.0f;
  bool  labelZTest = false;

  LambdaDict funcs;
  funcs["--obj"       ] = ArgFunc{1, [&](const char **args){ objPath = args[0]; }};
  funcs["--size"      ] = ArgFunc{2, [&](const char **args){ winW      = atoi(args[0]); winH = atoi(args[1]); }};
  funcs["--label-size"] = ArgFunc{1, [&](const char **args){ labelSize = atof(args[0]); }};
  funcs["--zoom"      ] = ArgFunc{1, [&](const char **args){ initZoom  = atof(args[0]); }};
  funcs["--fit"       ] = ArgFunc{0, [&](const char **args){ doAutoFit = true; }};
  funcs["--line-width"] = ArgFunc{1, [&](const char **args){ lineWidth = atof(args[0]); }};
  funcs["--label-ztest"] = ArgFunc{0, [&](const char **args){ labelZTest = true; }};

  process_args(argc, argv, funcs, false);

  // Handle negative window size (relative to display)
  if (winW < 0 || winH < 0) {
      if (SDL_Init(SDL_INIT_VIDEO) >= 0) {
          SDL_DisplayMode dm;
          if (SDL_GetCurrentDisplayMode(0, &dm) == 0) {
              if (winW < 0) winW = dm.w + winW;
              if (winH < 0) winH = dm.h + winH;
          }
          // We don't quit SDL here because AppSDL2OGL3 will call Init again, 
          // but SDL_Init is reference counted so it's fine.
      }
  }

  MeshViewerApp app(winW, winH);
  
  // Apply settings
  app.scene->renderer.setLineWidth(lineWidth);
  app.scene->text.charW = labelSize;
  app.scene->text.charH = labelSize;
  app.scene->labelZTest = labelZTest;
  app.screen->cam.zoom = initZoom;

  if (app.scene) {
    app.scene->loadOBJ(objPath.c_str());

    if (doAutoFit) {
        // Calculate bounding box
        Vec3f pmin = {1e30f, 1e30f, 1e30f};
        Vec3f pmax = {-1e30f, -1e30f, -1e30f};
        for (const auto& p : app.scene->mesh.points) {
            Vec3f pf = (Vec3f)p;
            pmin.x = fmin(pmin.x, pf.x); pmin.y = fmin(pmin.y, pf.y); pmin.z = fmin(pmin.z, pf.z);
            pmax.x = fmax(pmax.x, pf.x); pmax.y = fmax(pmax.y, pf.y); pmax.z = fmax(pmax.z, pf.z);
        }
        
        Vec3f center = (pmin + pmax) * 0.5f;
        Vec3f dim = pmax - pmin;
        float radius = dim.norm() * 0.5f;

        // Center model
        app.scene->modelPos = center;
        
        // Adjust camera
        // Heuristic: Place camera at distance proportional to radius
        // and set zoom proportional to radius
        app.screen->camDist = -radius * 3.0f;
        app.screen->cam.zoom = radius; 
        
        printf("Auto-fit: Center(%g %g %g) Radius %g -> CamDist %g Zoom %g\n", 
            center.x, center.y, center.z, radius, app.screen->camDist, app.screen->cam.zoom);
    }
  }

  app.loop(1000000);
  return 0;
}
