
#include "GUI.h" // THE HEADER

//Vec2i GUI_mouse_old_pos = (Vec2i){0,0};
//GUIAbstractPanel* GUI_mouse_panel = 0;
int GUI_fontTex = 0;

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
        //Draw::drawText( inputText.c_str(), fontTex, textSize, 0, 0 );
        Draw::drawText( inputText.c_str(), fontTex, textSize, 0 );
        Draw3D::drawLine( (Vec3f){curPos*textSize,0.0,0.0}, (Vec3f){curPos*textSize,textSize*2,0.0} );
    glPopMatrix();
}

void GUITextInput::viewHUD( const Vec2i& pos, int fontTex ){
    //printf("DEBUG GUITextInput::viewHUD \n");
    //Draw3D::drawText( inputText.c_str(), pos, fontTex, textSize, 0, 0 );
    //glDisable    ( GL_LIGHTING   );
    //glDisable    ( GL_DEPTH_TEST );
    //glShadeModel ( GL_FLAT       );
    glPushMatrix();
        glTranslatef( pos.x, pos.y, 0.0 );
        //Draw::billboardCam();
        //Draw::drawText( inputText.c_str(), fontTex, textSize, 0, 0 );
        Draw::drawText( inputText.c_str(), fontTex, fontSizeDef, 0 );
        Draw3D::drawLine( (Vec3f){curPos*fontSizeDef,0.0,0.0}, (Vec3f){curPos*fontSizeDef,fontSizeDef*2,0.0} );
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

void GUIAbstractPanel::view  ( ){
    // Draw2D::drawPointCross({xmin,ymax},5);
    glCallList( gllist );
};

void GUIAbstractPanel::moveBy(int dx, int dy){
    xmin += dx; ymin += dy;
    xmax += dx; ymax += dy;
    redraw=true; //tryRender();
};

void GUIAbstractPanel::moveTo(int x, int y){
    moveBy(x-xmin,y-ymin);
    //xmax = x+(xmax-xmin);
    //ymax = y+(ymax-ymin);
    //xmin = x;
    //ymin = y;
    //redraw=true; tryRender();
};

void GUIAbstractPanel::render(){
    glDisable   ( GL_LIGHTING    );
    glDisable   ( GL_DEPTH_TEST  );
    glShadeModel( GL_FLAT        );
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
    if(caption.length()>0) {
        Draw  ::setRGB( textColor );
        //int nchr = strlen(caption);
        //Draw2D::drawText( caption, nchr, {xmin, ymax-fontSizeDef*2},  0.0, GUI_fontTex, fontSizeDef );
        Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2},  0.0, GUI_fontTex, fontSizeDef );
    }
}

void GUIAbstractPanel::tryRender(){
    if(!redraw) return;
    if(gllist)glDeleteLists(gllist,1);
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
    render();
    glEndList();
    redraw=false;
};

void GUIAbstractPanel::initPanel( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ ){
    caption=caption_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_;
    redraw = true;
    //val_text = new char[nCharMax];
};

GUIAbstractPanel* GUIAbstractPanel::onMouse( int x, int y, const SDL_Event& event, GUI& gui ){ if( check(x,y) ){return this; }else{ return NULL; } };
void GUIAbstractPanel::onKeyDown( const SDL_Event& e, GUI& gui ){};
void GUIAbstractPanel::onText( const SDL_Event& e ){};

// ==============================
//       class GUIPanel
// ==============================

void GUIPanel::view ( ){
    // Draw2D::drawPointCross({xmin,ymax},5);
    glCallList( gllist );
    int xcur = xmin + curPos*fontSizeDef;
    Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+fontSizeDef*2} );
};

void GUIPanel::render(){
    glDisable( GL_LIGHTING );
    glDisable( GL_DEPTH_TEST);
    glShadeModel( GL_FLAT     );
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
    if(isSlider){
        Draw::setRGB(barColor);
        //Draw2D::drawRectangle ( xmin, ymax-2 xmin+val2x(value), ymax, true );
        //Draw2D::drawRectangle ( xmin, ymax-2*fontSizeDef, xmin+val2x(value), ymax, true );
        Draw2D::drawRectangle ( xmin, ymax-4*fontSizeDef, xmin+val2x(value), ymax-2*fontSizeDef, true );
    }
    Draw  ::setRGB( textColor );
    //int nch = strlen(caption);
    //Draw2D::drawText( caption, nch, {xmin, ymin+fontSizeDef*2,}, 0.0,  GUI_fontTex, fontSizeDef );
    Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2}, 0.0,  GUI_fontTex, fontSizeDef );
    int nch = inputText.length();
    if( nch > 0 ){
        Draw2D::drawText( inputText.c_str(), nch, {xmin, ymin}, 0.0, GUI_fontTex, fontSizeDef );
    }
}

/*
void GUIPanel::tryRender(){
    if(!redraw) return;
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
    render();
    glEndList();
    redraw=false;
};
*/

void GUIPanel::onKeyDown( const SDL_Event&  e ){
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

void GUIPanel::onText( const SDL_Event&  e ){
    if( SDL_GetModState() & KMOD_CTRL ) return;
    //char ch = e.text.text[0];
    //printf( "input event >>%s<<\n", e.text.text );
    //if((curPos+1)<nCharMax)curPos++;
    //val_text[curPos] = ch;
    //inputText.push_back(ch);
    inputText.insert(curPos,e.text.text); curPos++;
    redraw = true;
}

GUIAbstractPanel* GUIPanel::onMouse( int x, int y, const SDL_Event& event, GUI& gui ){
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
//     class  ScisorBox
// ==============================

//void MultiPanel::initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ ){
void ScisorBox::initScisor( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ ){
    caption=caption_;
    //xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_, ymax=ymax_; //fontTex=fontTex_;
    redraw = true;
};

void ScisorBox::apply( ){
    glEnable(GL_SCISSOR_TEST);
    glScissor(xmin,ymin,xmax-xmin,ymax-ymin);
}

void ScisorBox::render( ){
    glDisable   ( GL_LIGHTING    );
    glDisable   ( GL_DEPTH_TEST  );
    glShadeModel( GL_FLAT        );
    Draw  ::setRGB( textColor );
    Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, false );
    if(caption.length()>0) {
        //Draw  ::setRGB( textColor );
        //int nchr = strlen(caption);
        //Draw2D::drawText( caption, nchr, {xmin, ymax-fontSizeDef*2},  0.0, GUI_fontTex, fontSizeDef );
        Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2},  0.0, GUI_fontTex, fontSizeDef );
    }
};

/*
void ScisorBox::tryRender( ){
    if(!redraw) return;
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
    render();
    glEndList();
    redraw=false;
};
*/

GUIAbstractPanel* ScisorBox::onMouse( int x, int y, const SDL_Event&  event, GUI& gui ){
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        int ycut = ymax-fontSizeDef*2;
        //printf( "y %i ycut %i \n", y, ycut  );
        if(y>ycut){
            active = this;
            gui.dragged = this;
        }
    }
    return active;
};

// ==============================
//     class  MultiPanel
// ==============================

//void MultiPanel::initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ ){
void MultiPanel::initMulti( const std::string& caption_, int xmin_, int ymin_, int xmax_, int dy_, int nsubs_ ){
    //xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    caption =caption_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_, dy=dy_; //fontTex=fontTex_;
    nsubs=nsubs_;
    subs = new GUIPanel*[nsubs_];
    //int dy = (ymax-ymin-fontSizeDef)/nsubs;
    ymax=ymin + dy*nsubs + fontSizeDef*2;
    int yi = ymax-dy-fontSizeDef*2;
    for(int i=0; i<nsubs; i++){
        char buf[16];
        sprintf(buf,"val_%i",i);
        subs[i] = new GUIPanel( buf, xmin,yi,xmax,yi+dy, true, true );
        //subs[i]->initPanel(xmin,yi,xmax,yi+dy);
        //subs[i]->caption = new char[16];
        //sprintf(subs[i]->caption,"val_%i",i);
        //char buf[16];
        //sprintf(buf,"val_%i",i);
        //subs[i]->caption = buf;
        yi-=dy;
    }
    redraw = true;
};

void MultiPanel::moveBy(int dx, int dy){
    xmin+=dx; xmax+=dx;
    ymin+=dy; ymax+=dy;
    redraw=true; //tryRender();
    for(int i=0;i<nsubs;i++){
        subs[i]->moveBy(dx,dy);
    }
};

void MultiPanel::open(){
    opened = true;
    ymin = ymax - fontSizeDef*2 - dy*nsubs;
    redraw=true; //tryRender();
};

void MultiPanel::close(){
    opened = false;
    ymin = ymax - fontSizeDef*2;
    redraw=true; //tryRender();
};

void MultiPanel::toggleOpen(){
    if(opened){close();}else{open();}
    //printf( "opened %i \n", opened );
};

void MultiPanel::view( ){
    glCallList( gllist );
    if(opened){
        for(int i=0; i<nsubs; i++){
            subs[i]->draw();
        }
    }
};

/*
void MultiPanel::tryRender( ){
    GUIAbstractPanel::tryRender();
    for(int i=0; i<nsubs; i++){
        subs[i]->tryRender();
    }
};
*/

void MultiPanel::render( ){
    GUIAbstractPanel::render();
    for(int i=0; i<nsubs; i++){
        subs[i]->render();
    }
};

GUIAbstractPanel* MultiPanel::onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ){
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        active = this;
        if(opened){
            for(int i=0; i<nsubs; i++){
                active = subs[i]->onMouse ( x, y, event, gui );
                if(active) return active;
            }
        }
        if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
            if(event.button.button == SDL_BUTTON_LEFT){
                //GUI_mouse_old_pos;
                gui.dragged = this;
                //printf("clicked on MultiPanel Title \n");
                if(event.button.clicks > 1 ){ // double click
                    toggleOpen();
                }
            }
        }
    }
    return active;
};

// ==============================
//     class  DropDownList
// ==============================

void DropDownList::initList( const std::string& caption_, int xmin_, int ymin_, int xmax_, int nSlots_ ){
    caption=caption_;
    nSlots=nSlots_,xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymin+2*fontSizeDef*(nSlots+1);
    redraw = true;
};

//void DropDownList ::view ( ){
//    glCallList( gllist );
//};

DropDownList* DropDownList::addItem(const std::string& label){
    labels.push_back(label);
    return this;
};

void DropDownList::open(){
    bOpened = true;
    ymin = ymax - fontSizeDef*2*(nSlots+1);
    redraw=true; //tryRender();
};

void DropDownList::close(){
    bOpened = false;
    ymin = ymax - fontSizeDef*2;
    redraw=true; //tryRender();
};

void DropDownList ::render(){
    glDisable   ( GL_LIGHTING    );
    glDisable   ( GL_DEPTH_TEST  );
    glShadeModel( GL_FLAT        );
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
    Draw  ::setRGB( textColor );
    //bOpened = false;
    //bOpened = true;
    if(bOpened){
        Draw  ::setRGB( 0x00FF00 );
        Draw2D::drawRectangle ( xmin, ymax-(iSelected+2)*(fontSizeDef*2), xmax, ymax-(iSelected+1)*(fontSizeDef*2), true );
        Draw  ::setRGB( textColor );
        if(caption.length()>0) {
            Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2},  0.0, GUI_fontTex, fontSizeDef );
        }
        for(int i=0; i<nSlots; i++){
            int iItem = i+iItem0;
            if( iItem<labels.size() ){
                Draw2D::drawText( labels[iItem].c_str(), labels[iItem].length(), {xmin, ymax-(i+2)*2*fontSizeDef}, 0.0, GUI_fontTex, fontSizeDef );
                //sprintf( labels[],"val%i",i);
                    //sprintf(subs[i]->caption,"val%i",i);
            }
        }
    }else{
        //if( iSelected<nItems ){
        //int nch = strlen(caption);
        Draw2D::drawText( labels[iSelected].c_str(), labels[iSelected].length(), {xmin, ymin}, 0.0, GUI_fontTex, fontSizeDef );
        //}
    }
};

/*
void DropDownList ::tryRender( ){
    if(!redraw) return;
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
    render();
    glEndList();
    redraw=false;
};
*/

GUIAbstractPanel* DropDownList::onMouse ( int x, int y, const SDL_Event& event, GUI& gui ){
    if( check( x, y ) ){
        //if( event.type == SDL_MOUSEBUTTONUP ){
        if( event.type == SDL_MOUSEBUTTONDOWN ){
            if(event.button.button == SDL_BUTTON_LEFT){
                if(bOpened){
                    int i = ((ymax-y)/(2*fontSizeDef)) - 1;
                    if( (i>0)&&(i<nSlots) ){
                        i += iItem0;
                        i=_min(i,labels.size()-1);
                        i=_max(i,0);
                        iSelected = i;
                    }
                    close();
                }else{
                    open();
                }
                redraw = true;
            }
        }
        return this;
    }
    return 0;
};


// ==============================
//     class  TreeView
// ==============================


void TreeView::view( ){
    glCallList( gllist );
    Draw  ::setRGB( 0x00FF00 );
    Draw2D::drawRectangle ( xmin+1, ymax-(iSelected+2)*(fontSizeDef*2), xmax-1, ymax-(iSelected+1)*(fontSizeDef*2), false );
};

void TreeView::initTreeView( const std::string& caption_, int xmin_, int ymin_, int xmax_, int nSlots_ ){
    caption=caption_;
    nSlots=nSlots_,xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymin+2*fontSizeDef*(nSlots+1);
    redraw = true;
};

void TreeView::render( ){
    updateLines();
    GUIAbstractPanel::render();
    //glColor3f(0.0,0.0,1.0); Draw2D::drawPointCross({xmax,ymax},5);
    //glColor3f(1.0,0.0,0.0); Draw2D::drawPointCross({xmin,ymin},5);
    Draw  ::setRGB( textColor );
    int yoff = ymax - 4*fontSizeDef;
    glDisable(GL_DEPTH_TEST);
    for( TreeViewTree* tr : lines ){
        if(yoff<=ymin) break;
        std::string& str = tr->content.caption;
        Draw2D::drawText( str.c_str(), str.length(), {xmin+2*fontSizeDef*tr->content.level, yoff}, 0.0, GUI_fontTex, fontSizeDef );
        yoff-=2*fontSizeDef;
    }
};

void TreeView::updateLines( TreeViewTree& node, int level ){
    //printf( "%i: updateLines[%i] : '%s' \n", nitems, level, node.content.caption.c_str() );
    //nitems++;
    node.content.level = level;
    lines.push_back( &node );
    if(node.content.open){
        for(TreeViewTree* tr: node.branches ){
            updateLines( *tr, level+1 );
        }
    }
}

GUIAbstractPanel* TreeView::onMouse( int x, int y, const SDL_Event& event, GUI& gui ){
    if( check( x, y ) ){
        //if( event.type == SDL_MOUSEBUTTONUP ){
        if( event.type == SDL_MOUSEBUTTONDOWN ){
            if(event.button.button == SDL_BUTTON_LEFT){
                int i = ((ymax-y)/(2*fontSizeDef)) - 1;
                if( (i>0)&&(i<nSlots) ){
                    i += iItem0;
                    i=_min(i,lines.size()-1);
                    i=_max(i,0);
                    iSelected = i;
                }
                if(event.button.clicks > 1 ){ // double click
                    if(iSelected<lines.size()){
                        lines[iSelected]->content.open^=true;
                        redraw=true;
                    }
                }
            }
        }
        return this;
    }
    return 0;
};


// ==============================
//    class GUI
// ==============================

GUIAbstractPanel* GUI::addPanel( GUIAbstractPanel* panel ){ panels.push_back(panel); return panels.back(); }

void GUI::onEvent( int mouseX, int mouseY, const SDL_Event& event ){
    GUIAbstractPanel* active;
    switch( event.type ){
        case SDL_KEYDOWN:
            //if(focused){ focused->onKeyDown( event ); }else{ txt.onKeyDown(  event ); }; break;
            if(focused){ focused->onKeyDown( event, *this ); }
            break;
        case SDL_TEXTINPUT:
            //if(focused){ focused->onText   ( event ); }else{ txt.onText   ( event );  }; break;
            if(focused){ focused->onText   ( event ); }
            break;
        case SDL_MOUSEBUTTONDOWN:
            active = NULL; focused=NULL;
            for(GUIAbstractPanel* panel: panels){
                active =  panel->onMouse( mouseX, mouseY, event, *this );
                if(active)focused=active;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            if(event.button.button == SDL_BUTTON_LEFT){
                dragged = 0;
            }
            break;
        case SDL_MOUSEMOTION:
            SDL_MouseMotionEvent* event_ = (SDL_MouseMotionEvent*)&event;
            //if(GUI_mouse_panel) GUI_mouse_panel->moveTo( GUI_mouse_panel->xmin+event->xrel, GUI_mouse_panel->ymin+event->yrel );
            if(dragged){
                //printf(" GUI_globalEventHandler  SDL_MOUSEMOTION  %i %i \n", event_->xrel, -event_->yrel );
                dragged->moveBy( event_->xrel, -event_->yrel );
            }
            break;
    };
};

void GUI::draw(){
    //glLineWidth(1.0);
    //glLineWidth(0.5);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    for(GUIAbstractPanel* panel: panels){ panel->draw(); }
    if(focused){
        Draw::setRGB(focused->textColor);
        Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);
    }
}
