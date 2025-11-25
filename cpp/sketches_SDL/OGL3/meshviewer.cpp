#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <GL/glew.h>
#include <SDL2/SDL.h>

#include "Mat3.h"
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

  void loadOBJ(const char *path) {
    mesh.fromFileOBJ(path);
    mesh.polygonsToTriangles(false);
    mesh.tris2normals(true);
    mesh.findEdges();

    renderer.uploadMesh_d(mesh.points.size(), mesh.triangles.size(),
                          (int *)&mesh.triangles[0], (double *)&mesh.points[0],
                          (double *)&mesh.normals[0]);

    Vec2i *evts = mesh.exportEdgeVs();
    renderer.uploadLines_d(mesh.points.size(), mesh.edges.size(), (int *)evts,
                           (double *)&mesh.points[0]);
    delete[] evts;

    // Build labels: one label per vertex with its index as text
    text.clear();
    int n = mesh.points.size();
    for (int i = 0; i < n; i++) {
      char buf[16];
      sprintf(buf, "%d", i);
      Vec3f p = (Vec3f)mesh.points[i];
      p.z += 0.1f; // small offset towards camera
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

    GLuint dummyFont = makeDummyFontTex(128);
    text.init(dummyFont, 0.5f, 0.5f);
    // text.numGlyphs   = 128;
    // text.glyphOffset = 0;
    checkGLError("TextRendererOGL3::init");
  }

  virtual void draw(Camera &cam) override {
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    renderer.setModelPos(modelPos);
    renderer.draw(cam);

    text.setModelPos(modelPos);
    text.draw(cam);
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

  LambdaDict funcs;
  funcs["--obj"] = ArgFunc{1, [&](const char **args) { objPath = args[0]; }};
  process_args(argc, argv, funcs, false);

  int W = 800;
  int H = 800;
  MeshViewerApp app(W, H);
  if (app.scene) {
    app.scene->loadOBJ(objPath.c_str());
  }

  app.loop(1000000);
  return 0;
}
