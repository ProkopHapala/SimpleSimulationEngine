
#include "GLObject.h" // THE HEADER


void GLObject::init(){
    // vertexes to GPU

    for(int i=0; i<nbuffs; i++){
        /*
        if( buffs[i].cbuff == NULL ) continue;
        glGenBuffers( 1, &buffs[i].vbo );
        glBindBuffer( buffs[i].id, buffs[i].vbo );
        glBufferData( buffs[i].target, nVert*buffs[i].dim * sizeof(GLfloat), buffs[i].cbuff, buffs[i].usage );
        */
        if( index_cbuff  ) {
            //indexes.target=GL_ELEMENT_ARRAY_BUFFER;
            //indexes.dtype=GL_UNSIGNED_INT;
            //buffs[i].toGPU( nInd );
            glGenBuffers(1, &index_vbo);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo );
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, nInd * sizeof(unsigned int), index_cbuff, GL_STATIC_DRAW);
        }
        if( buffs[i].cbuff ) buffs[i].toGPU( nVert );
    }

    /*
    glGenBuffers(2, vbo); 						// Allocate and assign two Vertex Buffer Objects to our handle
    glBindBuffer(GL_ARRAY_BUFFER, vbo[0]); 		// Bind our first VBO as being the active buffer and storing vertex attributes (coordinates)
    glBufferData(GL_ARRAY_BUFFER, nVert*vboDim[0] * sizeof(GLfloat), vertexes, GL_STATIC_DRAW); 	// Copy the vertex data from diamond to our buffer  8 * sizeof(GLfloat) is the size of the diamond array, since it contains 8 GLfloat values
    // colors to GPU
    if(colors){
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]); 		                // Bind our second VBO as being the active buffer and storing vertex attributes (colors)
        glBufferData(GL_ARRAY_BUFFER, nVert*vboDim[1]*sizeof(GLfloat), colors, GL_STATIC_DRAW); // Copy the color data from colors to our buffer  12 * sizeof(GLfloat) is the size of the colors array, since it contains 12 GLfloat values
    }
    */
}

void GLObject::preDraw(){
    //printf(" preDraw( \n");
    for(int i=0; i<nbuffs; i++){ if (buffs[i].vbo ) buffs[i].activate(); }
}

void GLObject::afterDraw(){
    //printf(" afterDraw \n");
    for(int i=0; i<nbuffs; i++){ if(buffs[i].vbo)glDisableVertexAttribArray(buffs[i].id);  }
}

void GLObject::draw_instance(){
    //printf(" draw_instance\n");
    if(index_vbo){
        glBindBuffer  ( GL_ELEMENT_ARRAY_BUFFER, index_vbo );
        glDrawElements( draw_mode, nInd,  GL_UNSIGNED_INT, (void*)0 );
    }else{
        glDrawArrays( draw_mode, 0, nVert);
    }
}

void GLObject::draw_default(){
    //printf("==== draw_default \n");
    preDraw();
    draw_instance();
    afterDraw();
}

void GLObject::draw(){
    draw_default();
}

void GLObject::destroy(){
    for(int i=0; i<nbuffs; i++){ if(buffs[i].vbo) glDeleteBuffers(1, &buffs[i].vbo); }
}
