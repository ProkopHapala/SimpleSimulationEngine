
#include "GUI.h" // THE HEADER

// ==============================
//    class GUITextInput
// ==============================

void GUITextInput::applyVal( float f ){
    switch(num_op){
        case 0                : value =f; break;
        case SDLK_KP_PLUS     : value+=f; break;
        case SDLK_KP_MINUS    : value-=f; break;
        case SDLK_KP_MULTIPLY : value*=f; break;
        case SDLK_KP_DIVIDE   : value/=f; break;
    }
    if(checkRange){
        if      (value<vmin) {value=vmin;}
        else if (value>vmax) {value=vmax;}
    };
    num_op = 0;
}

void GUITextInput::view3D( const Vec3d& pos, int fontTex, float textSize ){
    //Draw3D::drawText( inputText.c_str(), pos, fontTex, textSize, 0, 0 );
    glDisable    ( GL_LIGHTING   );
    glDisable    ( GL_DEPTH_TEST );
    glShadeModel ( GL_FLAT       );
    glPushMatrix();
        glTranslatef( pos.x, pos.y, pos.z );
        Draw::billboardCam( );
        Draw::drawText( inputText.c_str(), fontTex, textSize, 0, 0 );
        Draw3D::drawLine( {curPos*textSize,0.0,0.0}, {curPos*textSize,textSize*2,0.0} );
    glPopMatrix();
}

void GUITextInput::onKeyDown( SDL_Event e ){
	// see https://wiki.libsdl.org/SDL_Keysym
    if ( SDL_GetModState() & KMOD_CTRL ){
        switch (e.key.keysym.sym ){
            case SDLK_v:
                inputText = SDL_GetClipboardText();        modified = true; break;
            case SDLK_c:
                SDL_SetClipboardText( inputText.c_str() ); modified = true; break;
        }
    }else{
        switch (e.key.keysym.sym ){
            case SDLK_BACKSPACE:
                if ( (inputText.length() > 0) && (curPos>0) ){ inputText.erase(curPos-1,1); curPos--; modified= true;} break;
            case SDLK_LEFT:
                if(curPos>0) curPos--; break;
            case SDLK_RIGHT:
                if(curPos<(inputText.length())) curPos++; break;
            case SDLK_RETURN:
            case SDLK_KP_ENTER:
                if( isNumber ){
                    try{
                        const char * str = inputText.c_str();
                        if(num_op!=0) str++; // skip numerical operation character
                        float f = std::stof( str );
                        applyVal( f );
                    }catch(std::exception const &exc){
                        printf("exception:%s\n", exc.what() );
                    };
                    //sprintf(inputText, "%f", value );
                    inputText = std::to_string(value);
                    num_op = 0;
                }
                entered=true;
                break;
            case SDLK_KP_PLUS:
            case SDLK_KP_MINUS:
            case SDLK_KP_DIVIDE:
            case SDLK_KP_MULTIPLY:
                if( isNumber ){
                    inputText = "";
                    curPos=0;
                    num_op = e.key.keysym.sym;
                    printf("num_op %i\n",num_op);
                }
                break;
        }
        printf("curPos : %i\n", curPos);
    }
};

void GUITextInput::onText( SDL_Event e ){
    if( SDL_GetModState() & KMOD_CTRL ) return;
    inputText.insert(curPos,e.text.text); curPos++;
    modified = true;
}

// ==============================
//    class GUIAbstractPanel
// ==============================

void GUIAbstractPanel::draw  ( ){
    glCallList( gllist );
};

void GUIAbstractPanel::tryRender(){
    if(!redraw) return;
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
        glDisable   ( GL_LIGHTING    );
        glDisable   ( GL_DEPTH_TEST  );
        glShadeModel( GL_FLAT        );
        Draw  ::setRGB( bgColor );
        Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
        printf("&caption %i\n", caption );
        if(caption) {
            Draw  ::setRGB( textColor );
            Draw2D::drawString ( caption,             xmin, ymax-12, 6, fontTex );
        }
    glEndList();
    redraw=false;
};

void GUIAbstractPanel::init( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_ ){
    xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    //val_text = new char[nCharMax];
};

// ==============================
//       class GUIPanel
// ==============================

void GUIPanel::draw  ( ){
    glCallList( gllist );
    int xcur = xmin + curPos*6;
    Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+12} );
};

void GUIPanel::tryRender(){
    //Draw3D::drawRect( xmin, ymin, xmax, ymax );
    if(!redraw) return;
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
        glDisable( GL_LIGHTING );
        glDisable( GL_DEPTH_TEST);
        glShadeModel( GL_FLAT     );
        //glColor3f   ( 0.8, 0.8, 0.8 );
        Draw  ::setRGB( bgColor );
        //printf("panel render %3.3f %3.3f %3.3f %3.3f \n", (float)xmin, (float)ymin, (float)xmax, (float)ymax );
        Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
        //float val_=(float)( (value-vmin)/(vmax-vmin) );
        if(isSlider){ Draw::setRGB(barColor); Draw2D::drawRectangle ( xmin, ymax-2, xmin+val2x(value), ymax, true ); }
        Draw  ::setRGB( textColor );
        Draw2D::drawString ( caption,             xmin, ymin+12, 6, fontTex );
        //Draw2D::drawString ( val_text, 0, curPos, xmin, ymin,    6, fontTex );
        int nch = inputText.length();
        if( nch > 0 ){
            Draw2D::drawString ( inputText.c_str(), 0, nch, xmin, ymin,    6, fontTex );
        }
    glEndList();
    redraw=false;
};

void GUIPanel::onKeyDown( SDL_Event e ){
    // see https://wiki.libsdl.org/SDL_Keysym
    if ( SDL_GetModState() & KMOD_CTRL ){
        switch (e.key.keysym.sym ){
            case SDLK_v:
                inputText = SDL_GetClipboardText(); redraw = true; break;
            case SDLK_c:
                SDL_SetClipboardText( inputText.c_str() ); redraw = true; break;
        }
    }else{
        switch (e.key.keysym.sym ){
            case SDLK_BACKSPACE:
                if ( (inputText.length() > 0) && (curPos>0) ){ inputText.erase(curPos-1,1); curPos--; redraw = true;} break;
            case SDLK_LEFT:
                if(curPos>0) curPos--; break;
            case SDLK_RIGHT:
                if(curPos<(inputText.length())) curPos++; break;
            case SDLK_RETURN:
            case SDLK_KP_ENTER:
                try{
                    float f = std::stof( inputText.c_str() );
                    value=f;
                }catch(std::exception const &exc){
                    printf("exception:%s\n", exc.what() );
                };
                executed = true;
                break;
        }
        printf("curPos : %i\n", curPos);
    }
};

void GUIPanel::onText( SDL_Event e ){
    if( SDL_GetModState() & KMOD_CTRL ) return;
    //char ch = e.text.text[0];
    //printf( "input event >>%s<<\n", e.text.text );
    //if((curPos+1)<nCharMax)curPos++;
    //val_text[curPos] = ch;
    //inputText.push_back(ch);
    inputText.insert(curPos,e.text.text); curPos++;
    redraw = true;
}

GUIAbstractPanel* GUIPanel::onMouse( int x, int y, SDL_Event event ){
    //printf( "panel.onMouse %i %i \n", x, y );
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        toRelative(x,y);
        //printf( "  panel.onMouse %i %i \n", x, y );
        if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
            active = this;
            if(isSlider && (event.button.button==SDL_BUTTON_RIGHT)){
                //value=( x*(vmax-vmin)/(xmax-xmin) ) + vmin;
                value=x2val(x);
                //sprintf(val_text, "%3.3f", value );
                inputText = std::to_string(value);
                redraw=true;
            }
            if(isButton && (event.button.button==SDL_BUTTON_LEFT ) ){
                executed=true;
                if (command!=NULL) command(value);
            }
        }
    }
    return active;
}

// ==============================
//     class  MultiPanel
// ==============================

void MultiPanel::initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ ){
    xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    nsubs=nsubs_;
    subs = new GUIPanel*[nsubs_];
    int dy = (ymax-ymin-6)/nsubs;
    int yi = 0;
    for(int i=0; i<nsubs; i++){
        subs[i] = new GUIPanel();
        subs[i]->init(xmin,yi,xmax,yi+dy, fontTex );
        subs[i]->caption = new char[16];
        sprintf(subs[i]->caption,"val%i",i);
        yi+=dy;
    }
};

void MultiPanel::draw  ( ){
    glCallList( gllist );
    for(int i=0; i<nsubs; i++){
        subs[i]->draw();
    }
};

void MultiPanel::tryRender( ){
    GUIAbstractPanel::tryRender();
    for(int i=0; i<nsubs; i++){
        subs[i]->tryRender();
    }
};

GUIAbstractPanel* MultiPanel::onMouse  ( int x, int y, SDL_Event event ){
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        for(int i=0; i<nsubs; i++){
            active = subs[i]->onMouse ( x, y, event );
            if(active) return active;
        }
    }
    return active;
};


