
#include "GLObject.h" // THE HEADER


void GLObject::init(){
    // vertexes to GPU

    for(int i=0; i<nbuffs; i++){
        if( buffs[i].cbuff == NULL ) continue;
        /*
        glGenBuffers( 1, &buffs[i].vbo );
        glBindBuffer( buffs[i].id, buffs[i].vbo );
        glBufferData( buffs[i].target, nVert*buffs[i].dim * sizeof(GLfloat), buffs[i].cbuff, buffs[i].usage );
        */
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

void GLObject::draw_default(){
    for(int i=0; i<nbuffs; i++){
        if (buffs[i].vbo ) buffs[i].activate();
        /*
        glEnableVertexAttribArray( buffs[i].id );
        glBindBuffer(GL_ARRAY_BUFFER, buffs[i].vbo );
        glVertexAttribPointer( buffs[i].id, buffs[i].dim, GL_FLOAT, GL_FALSE, 0, 0 );
        */
    }
    glDrawArrays( draw_mode, 0, nVert);
    for(int i=0; i<nbuffs; i++){ if(buffs[i].vbo)glDisableVertexAttribArray(buffs[i].id);  }

/*
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
    glVertexAttribPointer(0, vertDim, GL_FLOAT, GL_FALSE, 0, 0);
    if(vbo[1]){
        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        glVertexAttribPointer(1, clrDim, GL_FLOAT, GL_FALSE, 0, 0);
    }
    if(vbo[2]){
        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        glVertexAttribPointer(2, clrDim, GL_FLOAT, GL_FALSE, 0, 0);
    }
    glDrawArrays( draw_mode, 0, nVert);
    glDisableVertexAttribArray(0);
    if(vbo[1])glDisableVertexAttribArray(1);
    if(vbo[2])glDisableVertexAttribArray(2);
    //glDisableVertexAttribArray(2);
    */
}

void GLObject::draw(  ){
    draw_default();
}

void GLObject::destroy(){
    for(int i=0; i<nbuffs; i++){ if(buffs[i].vbo) glDeleteBuffers(1, &buffs[i].vbo); }
}
