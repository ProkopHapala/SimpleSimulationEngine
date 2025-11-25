#include "MeshRenderOGL3.h"

#include "IO_utils.h"

void MeshRenderOGL3::initDefaultShaders() {
  const char *namesShade[] = {"shade3D.vert", "shade3D.frag"};
  char **srcShade = fileGetSections("common_resources/shaders/Basic.glslf", 2,
                                    (char **)namesShade, (char *)"//>>");
  sh_solid.init_str(srcShade[0], srcShade[1], nullptr);
  sh_solid.getDefaultUniformLocation();

  sh_const.init_default();
  sh_const.getDefaultUniformLocation();

  modelPos.set(0.0f, 0.0f, 0.0f);
  modelMat.setOne();
}

void MeshRenderOGL3::setModelPos(const Vec3f &pos) { modelPos = pos; }

void MeshRenderOGL3::setModelMat(const Mat3f &mat) { modelMat = mat; }

void MeshRenderOGL3::uploadMesh_d(int nVerts, int nTris, const int *tris,
                                  const double *verts, const double *nors) {
  if (!mesh_tri)
    mesh_tri = new GLMesh();
  mesh_tri->init_d(nVerts, nTris * 3, (int *)tris, (double *)verts,
                   (double *)nors, nullptr, nullptr);
  mesh_tri->draw_mode = GL_TRIANGLES;

  // Create/Recreate VAO for this mesh configuration
  if (vaoTri)
    glDeleteVertexArrays(1, &vaoTri);
  glGenVertexArrays(1, &vaoTri);
  glBindVertexArray(vaoTri);

  // Position
  if (mesh_tri->buffs[0]) {
    glBindBuffer(GL_ARRAY_BUFFER, mesh_tri->buffs[0]);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
  }

  // Normal
  if (mesh_tri->buffs[1]) {
    glBindBuffer(GL_ARRAY_BUFFER, mesh_tri->buffs[1]);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
  }

  // Indices
  if (mesh_tri->inds) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri->inds);
  }

  glBindVertexArray(0);
}

void MeshRenderOGL3::uploadLines_d(int nVerts, int nEdges, const int *edges,
                                   const double *verts) {
  if (!mesh_lines)
    mesh_lines = new GLMesh();
  mesh_lines->init_d(nVerts, nEdges * 2, (int *)edges, (double *)verts, nullptr,
                     nullptr, nullptr);
  mesh_lines->draw_mode = GL_LINES;

  // Create/Recreate VAO for this mesh configuration
  if (vaoLine)
    glDeleteVertexArrays(1, &vaoLine);
  glGenVertexArrays(1, &vaoLine);
  glBindVertexArray(vaoLine);

  // Position
  if (mesh_lines->buffs[0]) {
    glBindBuffer(GL_ARRAY_BUFFER, mesh_lines->buffs[0]);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);
  }

  // Indices
  if (mesh_lines->inds) {
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_lines->inds);
  }

  glBindVertexArray(0);
}

void MeshRenderOGL3::draw(const Camera &cam) {
  // Build camera matrix from Camera parameters (similar to setCameraPersp in
  // SceneOGL3.h)
  Mat4f camMat, mRot, mPersp;
  if (cam.persp) {
    mPersp.setPerspective(cam.zoom / cam.aspect, cam.zoom, cam.zmin, cam.zmax);
    mRot.setOne();
    mRot.setRot(cam.rot);
    camMat.set_mmul_TN(mRot, mPersp);
  } else {
    camMat.setOrthographic(cam.zoom, cam.zoom * cam.aspect, cam.zmin, cam.zmax);
  }

  // Debug: print camera matrix
  // printf("MeshRenderOGL3::draw() camMat: \n"); camMat.print();

  if (mesh_tri && vaoTri) {
    sh_solid.use();
    sh_solid.set_modelPos((GLfloat *)&modelPos);
    sh_solid.set_modelMat((GLfloat *)&modelMat);
    sh_solid.set_camPos((GLfloat *)&cam.pos);
    sh_solid.set_camMat((GLfloat *)&camMat);

    // Lighting setup (hardcoded defaults for now)
    glUniform3fv(sh_solid.getUloc((char*)"light_pos"),     1, (const float[]){ 1.0f,  1.0f, -1.0f });
    glUniform3fv(sh_solid.getUloc((char*)"lightColor"),    1, (const float[]){ 1.0f,  0.9f,  0.8f });
    glUniform3fv(sh_solid.getUloc((char*)"diffuseColor"),  1, (const float[]){ 1.0f,  1.0f,  1.0f });
    glUniform3fv(sh_solid.getUloc((char*)"ambientColor"),  1, (const float[]){ 0.2f,  0.2f,  0.3f });
    glUniform3fv(sh_solid.getUloc((char*)"specularColor"), 1, (const float[]){ 0.0f,  0.0f,  0.0f });

    glBindVertexArray(vaoTri);
    glDrawElements(GL_TRIANGLES, mesh_tri->nInds, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
  }
  if (mesh_lines && vaoLine) {
    sh_const.use();
    sh_const.set_modelPos((GLfloat *)&modelPos);
    sh_const.set_modelMat((GLfloat *)&modelMat);
    sh_const.set_camPos((GLfloat *)&cam.pos);
    sh_const.set_camMat((GLfloat *)&camMat);
    const float colRed[4] = {1.0f, 0.0f, 0.0f, 1.0f};
    const float colGreen[4] = {0.0f, 1.0f, 0.0f, 1.0f};
    sh_const.set_baseColor(colGreen);
    glLineWidth(lineWidth);

    glBindVertexArray(vaoLine);
    glDrawElements(GL_LINES, mesh_lines->nInds, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
  }
}
