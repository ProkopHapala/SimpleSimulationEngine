
#define X .525731112119133606
#define Z .850650808352039932

void drawtriangle(float *v1, float *v2, float *v3){
  glBegin(GL_TRIANGLES);
    glNormal3fv(v1);    glVertex3fv(v1);
    glNormal3fv(v2);    glVertex3fv(v2);
    glNormal3fv(v3);    glVertex3fv(v3);
    glFrontFace(GL_CW);
  glEnd();
}

void normalize(float v[3]){
  float d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d == 0.0){
    printf("zero length vector");
    return;
  }
  v[0] /= d; v[1] /= d; v[2] /= d;
}

void subdivide(float *v1, float *v2, float *v3, long depth){
  float v12[3], v23[3], v31[3];
  if(depth <= 0){
    drawtriangle(v1, v2, v3);
    return;
  }
  for(int i=0; i<3; i++){
    v12[i] = (v1[i]+v2[i])*0.5;
    v23[i] = (v2[i]+v3[i])*0.5;
    v31[i] = (v3[i]+v1[i])*0.5;
  }
  normalize(v12);
  normalize(v23);
  normalize(v31);
  subdivide( v1, v12, v31, depth-1);
  subdivide( v2, v23, v12, depth-1);
  subdivide( v3, v31, v23, depth-1);
  subdivide(v12, v23, v31, depth-1);
}

void normcrossprod(float v1[3], float v2[3], float out[3]){
  out[0] = v1[1]*v2[2] - v1[2]*v2[1];
  out[1] = v1[2]*v2[0] - v1[0]*v2[2];
  out[2] = v1[0]*v2[1] - v1[1]*v2[0];
  normalize(out);
}

// ====== sisplay func

void drawIcosaSphere( int level ){
  static float vdata[12][3] = {
    {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
    {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
    {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}  
  };
  static int tindices[20][3] = {
    {1,4,0},  {4,9,0},  {4,5,9},  {8,5,4},  {1,8,4},
    {1,10,8}, {10,3,8}, {8,3,5},  {3,2,5},  {3,7,2},
    {3,10,7}, {10,6,7}, {6,11,7}, {6,0,11}, {6,1,0}, 
    {10,1,6}, {11,0,9}, {2,11,9}, {5,2,9},  {11,2,7},   
  };
  float d1[3], d2[3], norm[3];
  glBegin(GL_TRIANGLES);
  for(int i=0; i<20; i++){
    for( int j=0; j<3; j++){
      d1[j] = vdata[tindices[i][0]][j] - vdata[tindices[i][1]][j];
      d2[j] = vdata[tindices[i][1]][j] - vdata[tindices[i][2]][j];
    }
    normcrossprod(d1, d2, norm);
    glNormal3fv(norm);
    subdivide(&vdata[tindices[i][0]][0], &vdata[tindices[i][1]][0], &vdata[tindices[i][2]][0],  level  );
  }
  glEnd();

}



