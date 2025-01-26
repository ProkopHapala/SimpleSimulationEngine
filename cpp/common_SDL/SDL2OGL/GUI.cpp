
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
        Draw3D::drawLine( Vec3f{curPos*textSize,0.0,0.0}, Vec3f{curPos*textSize,textSize*2,0.0} );
    glPopMatrix();
}

void GUITextInput::viewHUD( const Vec2i& pos, int fontTex, bool bBack ){
    glPushMatrix();
        glTranslatef( pos.x, pos.y, 0.0 );
        //Draw::billboardCam();
        //Draw::drawText( inputText.c_str(), fontTex, textSize, 0, 0 );
        int nl = inputText.size();
        if(bBack)Draw2D::drawRectangle( (Vec2f){pos.x,pos.y}, (Vec2f){pos.x+nl*fontSizeDef, pos.y+fontSizeDef*2}, true );
        Draw::drawText( inputText.c_str(), fontTex, fontSizeDef, 0 );
        Draw3D::drawLine( Vec3f{curPos*fontSizeDef,0.0,0.0}, Vec3f{curPos*fontSizeDef,fontSizeDef*2,0.0} );
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
                    //printf("num_op %i\n",num_op);
                }
                break;
        }
        //printf("curPos : %i\n", curPos);
    }
}

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
}

void GUIAbstractPanel::moveBy(int dx, int dy){
    xmin += dx; ymin += dy;
    xmax += dx; ymax += dy;
    redraw=true; //tryRender();
}

void GUIAbstractPanel::moveTo(int x, int y){
    moveBy(x-xmin,y-ymin);
    //xmax = x+(xmax-xmin);
    //ymax = y+(ymax-ymin);
    //xmin = x;
    //ymin = y;
    //redraw=true; tryRender();
}

void GUIAbstractPanel::render(){
    //printf( "GUIAbstractPanel::render() p0(%i,%i) p2(%i,%i) \n", xmin, ymin, xmax, ymax );
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
    redraw=false;
}

void GUIAbstractPanel::tryRender(){
    if(!redraw) return;
    if(gllist)glDeleteLists(gllist,1);
    gllist=glGenLists(1);
    glNewList( gllist, GL_COMPILE );
    render();
    glEndList();
    redraw=false;
}

void GUIAbstractPanel::initPanel( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ ){
    //printf( "GUIAbstractPanel::initPanel(%s,pmin(%i,%i),pmax(%i,%i)) \n", caption_.c_str(), xmin_, ymin_, xmax_, ymax_ );
    caption=caption_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_;
    redraw = true;
    //val_text = new char[nCharMax];
}

GUIAbstractPanel* GUIAbstractPanel::onMouse( int x, int y, const SDL_Event& event, GUI& gui ){ if( check(x,y) ){return this; }else{ return NULL; } };
void GUIAbstractPanel::onKeyDown( const SDL_Event& e, GUI& gui ){ printf("GUIAbstractPanel::onKeyDown() not implemented\n"); };
void GUIAbstractPanel::onText   ( const SDL_Event& e, GUI& gui ){ printf("GUIAbstractPanel::onText() not implemented\n"); };

// ==============================
//       class GUIPanel
// ==============================

bool GUIPanel::checkRange(bool bExit, bool bWarn){
    bool ret = false;
    if( vmin>vmax                               ){ if(bWarn){ printf("WARRNING GUIPanel(%s) vmin(%g)>vmax(%g)\n", caption.c_str(), vmin, vmax );  }; ret=true; }
    if((vmax-vmin)<1e-8*(fabs(vmax)+fabs(vmin)) ){ if(bWarn){ printf("WARRNING GUIPanel(%s) (vmax(%g)-vmin(%g))=%g is numerically unstable \n", caption.c_str(), vmin, vmax, vmax-vmin );  }; ret=true; }
    if(ret && bExit                             ){            printf("ERROR in GUIPanel(%s)::checkRange()=>exit()", caption.c_str() ); exit(0);}
    return ret;
}

bool GUIPanel::checkValue(bool bExit, bool bWarn){
    bool ret = false;
    if( (value>vmax)||(value<vmin) ){ if(bWarn){ printf("WARRNING GUIPanel(%s) value(%g) out of range vmin(%g)..vmax(%g)\n", caption.c_str(), value, vmin, vmax );  }; ret=true; }
    if(ret && bExit){ printf( "ERROR in GUIPanel(%s)::checkValue()=>exit()", caption.c_str() ); exit(0);}
    return ret;
}

void GUIPanel::view ( ){
    //tryRender();
    //Draw2D::drawPointCross({xmin,ymax},5);
    glCallList( gllist );
    int nch0 = caption.length();
    int xcur = xmin + (nch0+curPos)*fontSizeDef;
    Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+fontSizeDef*2} );
}

void GUIPanel::render(){
    //printf( "GUIPanel(%s)::render() this(%li) p0(%i,%i) p2(%i,%i) isSlider(%i) isButton(%i) \n", caption.c_str(), (long)this, ymin, xmax, ymax, isSlider, isButton );
    if(isInt){ value=getIntVal(); }
    glDisable( GL_LIGHTING   );
    glDisable( GL_DEPTH_TEST );
    glShadeModel( GL_FLAT    );
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );

    // Border ?
    //Draw  ::setRGB( textColor ); Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, false );

    if(isSlider){
        Draw::setRGB(barColor);
        float dx = val2x(value);
        if( (dx<0)||(dx>(xmax-xmin)) ){ printf( "GUIPanel()::render() va(%f) out of range(0,%g) \n", caption.c_str(), dx, (xmax-xmin) ); }
        Draw2D::drawRectangle ( xmin, ymin, xmin+dx, ymax, true );
        //Draw2D::drawRectangle ( xmin, ymax-2*fontSizeDef, xmin+val2x(value), ymax, true );
        //Draw2D::drawRectangle ( xmin, ymin, xmin+val2x(value), ymax-2*fontSizeDef, true );
    }
    Draw  ::setRGB( textColor );
    int nch0 = caption.length();
    //Draw2D::drawText( caption, nch, {xmin, ymin+fontSizeDef*2,}, 0.0,  GUI_fontTex, fontSizeDef );
    Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2}, 0.0,  GUI_fontTex, fontSizeDef );
    if(viewVal){ val2text(); }
    int nch = inputText.length();
    if( nch > 0 ){
        //Draw  ::setRGB( 0xFFFFFFFF );
        //Draw2D::drawRectangle( xmin+nch0*fontSizeDef, ymax-2*fontSizeDef, xmax, ymax, true );
        Draw  ::setRGB( textColor );
        Draw2D::drawText( inputText.c_str(), nch, {xmin+fontSizeDef*nch0, ymin}, 0.0, GUI_fontTex, fontSizeDef );
    }
    redraw=false;
}

void GUIPanel::onKeyDown( const SDL_Event&  e, GUI& gui  ){
    bool doIt = false;
    //printf( "GUIPanel(%s)::onKeyDown() key=%i \n", caption.c_str(), e.key.keysym.sym );
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
            case SDLK_DELETE:
                if ( (inputText.length() > 0) && (curPos<inputText.length()) ){ inputText.erase(curPos,1); redraw = true;} break;
            case SDLK_LEFT:
                if(curPos>0) curPos--; break;
            case SDLK_RIGHT:
                if(curPos<(inputText.length())) curPos++; break;
            case SDLK_RETURN:
            case SDLK_KP_ENTER: doIt = true; [[fallthrough]];
            case SDLK_TAB:
                try{
                    float f;
                    if(isInt){ f=getIntVal(); }else{ std::stof( inputText.c_str() ); }
                    value=f;
                    redraw=true;
                }catch(std::exception const &exc){  
                    //printf("exception:%s\n", exc.what() ); 
                    printf("GUIPanel(%s)::onKeyDown() problem convert inputText(%s) to value(%g) | exception:%s\n", caption.c_str(), inputText.c_str(), value, exc.what() );
                };
                if(doIt) command(this);
                executed = true;
                break;
        }
        //printf("curPos : %i\n", curPos);
    }
}

void GUIPanel::onText( const SDL_Event&  e, GUI& gui ){
    //printf( "GUIPanel(%s)::onText() text=%s \n", caption.c_str(), e.text.text );
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
        //printf( "  panel.onMouse %i %i isSlider=%i \n", x, y, isSlider );
        if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
            active = this;
            if(isSlider && (event.button.button==SDL_BUTTON_RIGHT)){ // onSliderChange 
                //printf( "  panel.onMouse SDL_BUTTON_RIGHT  \n" );
                //value=( x*(vmax-vmin)/(xmax-xmin) ) + vmin;
                ivalchanged=1;
                value=x2val(x);
                if(isInt){ value=getIntVal(); }
                if(bCmdOnSlider) command(this);
                //sprintf(val_text, "%3.3f", value );
                //inputText = std::to_string(value);
                redraw=true;
            }
            if(isButton && (event.button.button==SDL_BUTTON_LEFT ) ){
                executed=true;
                //if (command!=NULL) command(value,caller);
                if (command) command(this);
            }else{
                SDL_StartTextInput();
            }
        }
    }
    return active;
}

// ==============================
//     class  MultiPanel
// ==============================

//void MultiPanel::initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ ){
void MultiPanel::initMulti( const std::string& caption_, int xmin_, int ymin_, int xmax_, int dy_, int nsubs_, bool isSlider, bool isButton, bool isInt, bool viewVal, bool bCmdOnSlider ){
    //printf( "MultiPanel::initMulti(%s,nsubs=%i,dy=%i,dx=%i) pmin(%i,%i)  isSlider=%i isButton=%i isInt=%i viewVal=%i bCmdOnSlider=%i\n", caption_.c_str(), nsubs_, dy_, xmax_-xmin_, xmin_, ymin_, isSlider, isButton, isInt, viewVal, bCmdOnSlider );
    //xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    caption =caption_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_, dy=dy_; //fontTex=fontTex_;
    //nsubs=nsubs_;
    //subs = new GUIPanel*[nsubs_];
    bool balloc = nsubs_>0; 
    if(balloc){nsubs=nsubs_; subs.resize(nsubs); }else{nsubs=-nsubs_; subs.reserve(nsubs); }
    //int dy = (ymax-ymin-fontSizeDef)/nsubs;
    ymax=ymin + dy*nsubs + fontSizeDef*2;
    int yi = ymax-dy-fontSizeDef*2;
    if(balloc)
    for(int i=0; i<nsubs; i++){
        char buf[16];
        sprintf(buf,"val_%i",i);
        //printf( "MultiPanel::initMulti(%i) p0(%i,%i) p1(%i,%i) \n", i, xmin,yi,xmax,yi+dy );
        subs[i] = new GUIPanel( buf, xmin,yi,xmax,yi+dy, isSlider, isButton,isInt, viewVal, bCmdOnSlider );
        yi-=dy;
    }
    redraw = true;
}

void MultiPanel::moveBy(int dx, int dy){
    xmin+=dx; xmax+=dx;
    ymin+=dy; ymax+=dy;
    redraw=true; //tryRender();
    for(int i=0;i<nsubs;i++){
        subs[i]->moveBy(dx,dy);
    }
}

void MultiPanel::open(){
    opened = true;
    ymin = ymax - fontSizeDef*2 - dy*nsubs;
    redraw=true; //tryRender();
}

void MultiPanel::close(){
    opened = false;
    ymin = ymax - fontSizeDef*2;
    redraw=true; //tryRender();
}

void MultiPanel::toggleOpen(){
    if(opened){close();}else{open();}
    //printf( "opened %i \n", opened );
}

void MultiPanel::view( ){
    if( visible==false ) return;
    //printf( "MultiPanel::view() opened %i \n", opened );
    glCallList( gllist );
    // --- NOTE: we do not need to call view() for subs, because they are already baked into gllist ( see MultiPanel::render() )
    if(opened){ for(int i=0; i<nsubs; i++){ redraw |= subs[i]->redraw;} }
    // if(opened){
    //     for(int i=0; i<nsubs; i++){
    //         subs[i]->tryRender();
    //         //subs[i]->view();
    //     }
    // }
    //printf( "MultiPanel::view() END \n" );
}

void MultiPanel::render( ){
    if( visible==false ) return;
    //printf( "MultiPanel::render() opened=%i \n", opened );
    nsubs = subs.size();
    GUIAbstractPanel::render();
    if(opened){
        for(int i=0; i<nsubs; i++){
            subs[i]->render();
        }
    }
    //printf( "MultiPanel::render() END\n" );
}

GUIAbstractPanel* MultiPanel::onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ){
    if( visible==false ) return 0;
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        active = this;
        if(opened){
            nsubs = subs.size();
            for(int i=0; i<nsubs; i++){
                active = subs[i]->onMouse ( x, y, event, gui );
                if(subs[i]->redraw) redraw = true;
                if(active){
                    if(hideOnCommand){ 
                        visible = false; redraw = true; 
                        return this;
                    }; 
                    return active;
                }
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
}

// ==============================
//     class  CheckBoxList
// ==============================

void CheckBoxList::initCheckBoxList( int xmin_, int ymin_, int xmax_, int dy_){
    xmin=xmin_,ymin=ymin_,xmax=xmax_, dy=dy_; //fontTex=fontTex_;
    //ymax=ymin + dy*boxes.s + fontSizeDef*2;
    //int yi = ymax-dy-fontSizeDef*2;
    redraw = true;
}

/*
void CheckBoxList::moveBy(int dx, int dy){
    xmin+=dx; xmax+=dx;
    ymin+=dy; ymax+=dy;
    redraw=true; //tryRender();
    for(int i=0;i<nsubs;i++){
        subs[i]->moveBy(dx,dy);
    }
};
*/

void CheckBoxList::view( ){
    glCallList( gllist );
}

void CheckBoxList::update(){
    ymax=ymin + dy*boxes.size() + fontSizeDef*2;
    syncRead();
}

void CheckBoxList::render( ){
    glDisable( GL_LIGHTING );
    glDisable( GL_DEPTH_TEST);
    glShadeModel( GL_FLAT );
    update();

    int y0 = ymin+boxes.size()*dy;
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( xmin, y0, xmax, y0+dy, true );
    Draw  ::setRGB( textColor );
    Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2}, 0.0,  GUI_fontTex, fontSizeDef );
    //Draw2D::drawText( caption.c_str(), 0, {xmin, y0}, 0.0, GUI_fontTex, fontSizeDef );
    for(int i=0; i<boxes.size(); i++){
        const CheckBox& box = boxes[i];
        if(box.val){ Draw::setRGB(checkColor); }else{ Draw::setRGB(bgColor); }
        int y=ymin+i*dy;
        Draw2D::drawRectangle ( xmin, y, xmax, y+dy, true );
        Draw  ::setRGB( textColor );
        Draw2D::drawText( box.label.c_str(), 0, {xmin, y}, 0.0, GUI_fontTex, fontSizeDef );
    }
}

GUIAbstractPanel* CheckBoxList::onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ){
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        active = this;
        if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
            if(event.button.button == SDL_BUTTON_LEFT){
                gui.dragged = this;
                int ibox = (y-ymin)/dy;
                if( ibox<boxes.size()){
                    boxes[ibox].flip();
                    redraw=true;
                    ivalchanged=ibox;
                    return active;
                }
                //if(event.button.clicks > 1 ){ toggleOpen();}
            }
        }
    }
    ivalchanged=-1;
    return active;
}

// ==============================
//     class  ScisorBox
// ==============================

void ScisorBox::initScisor( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ ){
    caption=caption_;
    //xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
    xmin=xmin_,ymin=ymin_,xmax=xmax_, ymax=ymax_; //fontTex=fontTex_;
    redraw = true;
}

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
}

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

GUIAbstractPanel* ScisorBox::onMouse( int x, int y, const SDL_Event&  event, GUI& gui){
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
}

// ==============================
//     class  CommandList
// ==============================

void CommandList::initCommandList( int xmin_, int ymin_, int xmax_, int dy_){
    xmin=xmin_,ymin=ymin_,xmax=xmax_, dy=dy_; //fontTex=fontTex_;
    redraw = true;
}

void CommandList::view( ){
    glCallList( gllist );
}

void CommandList::update(){
    ymax=ymin + dy*commands->commands.size() + fontSizeDef*2;
}

bool CommandList::getKeyb(int key){
    bool doit = (icmdbind !=-1);
    if( doit ){ commands->rebind( icmdbind, key ); redraw = true; }
    icmdbind = -1;
    return doit;
}

void CommandList::render( ){
    glDisable( GL_LIGHTING );
    glDisable( GL_DEPTH_TEST);
    glShadeModel( GL_FLAT );
    update();
    auto& cmds = commands->commands;
    char stmp[256];
    for(int i=0; i<cmds.size(); i++){
        const Command& cmd = cmds[i];
        //if(box.val){ Draw::setRGB(checkColor); }else{ Draw::setRGB(bgColor); }
        int y=ymin+i*dy;
        if(i==icmdbind){ Draw::setRGB( modColor ); }else{ Draw::setRGB( bgColor ); };
        Draw2D::drawRectangle ( xmin, y, xmax, y+dy, true );
        Draw  ::setRGB( textColor );
        if( (cmd.key>=32)&&(cmd.key<128) ){ sprintf( stmp, " '%c' %s" , (char)cmd.key, cmd.name.c_str() ); }
        else                              { sprintf( stmp, "#%03i %s",       cmd.key, cmd.name.c_str() ); }
        Draw2D::drawText( stmp, 0, {xmin, y}, 0.0, GUI_fontTex, fontSizeDef );
    }
}

GUIAbstractPanel* CommandList::onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ){
    GUIAbstractPanel* active = NULL;
    if( check( x, y ) ){
        active = this;
        if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
            if(event.button.button == SDL_BUTTON_LEFT){
                //gui.dragged = this;
                int ibox = (y-ymin)/dy;
                auto& cmds = commands->commands;
                if( ibox<cmds.size()){
                    //if(commandDispatch) (*commandDispatch)( ibox );
                    if(bDispatch) commandDispatch( ibox );
                }
                //if(event.button.clicks > 1 ){ toggleOpen();}
            }else
            if(event.button.button == SDL_BUTTON_RIGHT){
                //gui.dragged = this;
                int ibox = (y-ymin)/dy;
                auto& cmds = commands->commands;
                if( ibox<cmds.size()){
                    icmdbind = ibox;
                    redraw = true;
                    // ToDo : modify keys
                    // Proble:  we need to capture next key
                }
                //if(event.button.clicks > 1 ){ toggleOpen();}
            }
        }
    }
    return active;
}



// ==============================
//     class  DropDownList
// ==============================

void DropDownList::initList( const std::string& caption_, int xmin_, int ymin_, int xmax_, int nSlots_ ){
    caption=caption_;
    nSlots=nSlots_,xmin=xmin_,ymin=ymin_,xmax=xmax_;//ymax=ymin+2*fontSizeDef*(nSlots+1);
    ymax=ymin+2*fontSizeDef;
    redraw = true;
}

int DropDownList::selectedToStr(char* str){
    //printf( "DropDownList::selectedToStr  %i '%s' \n", iSelected, labels[iSelected].c_str()  );
    return sprintf(str,"%s", labels[iSelected].c_str() );
}

//void DropDownList ::view ( ){
//    glCallList( gllist );
//}

DropDownList* DropDownList::addItem(const std::string& label){
    labels.push_back(label);
    return this;
}

void DropDownList::open(){
    bOpened = true;
    ymin = ymax - fontSizeDef*2*(nSlots+1);
    redraw=true; //tryRender();
}

void DropDownList::close(){
    bOpened = false;
    ymin = ymax - fontSizeDef*2;
    redraw=true; //tryRender();
}

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
        int icur = iSelected-iItem0;
        if((icur>=0)&&(icur<nSlots)) Draw2D::drawRectangle ( xmin, ymax-(icur+2)*(fontSizeDef*2), xmax, ymax-(icur+1)*(fontSizeDef*2), true );
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
}

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
    //printf( "DropDownList::onMouse x,y, %i(%i..%i) %i(%i..%i) \n", x,xmin,xmax,  y,ymin,ymax );
    //if( event.type == SDL_MOUSEWHEEL ){ printf( " SDL_MOUSEWHEEL !!! \n" ); };
    if( check( x, y ) ){
        //printf( "DropDownList::onMouse inside \n" );
        //if( event.type == SDL_MOUSEBUTTONUP ){
        if( event.type == SDL_MOUSEBUTTONDOWN ){
            if(event.button.button == SDL_BUTTON_LEFT){
                if(bOpened){
                    int i = ((ymax-y)/(2*fontSizeDef))-1;
                    //printf( "i %i \n", i );
                    if( (i>=0)&&(i<nSlots) ){
                        i += iItem0;
                        i=_min(i,(int)labels.size()-1);
                        i=_max(i,0);
                        iSelected = i;
                        //if(onSelect)onSelect->GUIcallback(this);
                        if(command)command(this);
                        //printf( "DropDownList::onMouse() iSelected %i  iItem0 %i  labels.size() %i  \n", iSelected, iItem0, labels.size() );
                    }
                    close();
                }else{
                    open();
                }
                redraw = true;
            }
        }else if( event.type == SDL_MOUSEWHEEL ){
            //printf( " SDL_MOUSEWHEEL \n" );
            if     (event.wheel.y < 0){ iItem0 = _min( iItem0+1, (int)labels.size()-nSlots ); }
            else if(event.wheel.y > 0){ iItem0 = _max( iItem0-1, 0                    ); }
            redraw = true;
        }
        return this;
    }
    return 0;
}


// ==============================
//     class  TreeView
// ==============================


void TreeView::view( ){
    glCallList( gllist );
    Draw  ::setRGB( 0x00FF00 );
    Draw2D::drawRectangle ( xmin+1, ymax-(iSelected+2)*(fontSizeDef*2), xmax-1, ymax-(iSelected+1)*(fontSizeDef*2), false );
}

void TreeView::initTreeView( const std::string& caption_, int xmin_, int ymin_, int xmax_, int nSlots_ ){
    caption=caption_;
    nSlots=nSlots_,xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymin+2*fontSizeDef*(nSlots+1);
    redraw = true;
}

void TreeView::render( ){
    //if(content.open){
        updateLines();
        ymax=ymin+2*fontSizeDef*(nSlots+1);
    //}else{
    //    ymax=ymin+2*fontSizeDef;
    //}
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
}

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
                    i=_min(i,(int)lines.size()-1);
                    i=_max(i,0);
                    iSelected = i;
                }
                if(event.button.clicks > 1 ){ // double click
                    if(iSelected<lines.size()){
                        lines[iSelected]->content.open^=true;
                        redraw=true;
                    }
                }
                gui.dragged = this;
            }
        }
        return this;
    }
    return 0;
}


// ==============================
//     class  TableView
// ==============================

void TableView::initTableView( Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_ ){
    table = table_;
    caption=caption_;
    xmin=xmin_,ymin=ymin_,
    //xmax=xmax_,ymax=ymax_;
    i0=i0_; j0=j0_; imax=imax_; jmax=jmax_;
    int nch    = 8;
    int nchpix = nch*fontSizeDef;
    nchs.resize( table->columns.size() );
    xs  .resize( table->columns.size()+1 );
    int x=0;
    xs[0]=0;
    for(int i=0; i<nchs.size(); i++){ nchs[i]=nch; x+=nchpix; xs[i+1]=x; }
    xmax = xmin + nchpix       *(jmax-j0);
    ymax = ymin + fontSizeDef*2*(imax-i0);
    //printf( " i (%i,%i) j (%i,%i) \n", i0, imax,    j0, jmax   );
    //printf( " x (%i,%i) y (%i,%i) \n", xmin, xmax,  ymin, ymax );
    redraw = true;
}

void TableView::render(){
    glDisable( GL_LIGHTING );
    glDisable( GL_DEPTH_TEST);
    glShadeModel( GL_FLAT );
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle( xmin, ymin, xmax, ymax, true );
    //int ncol = table->columns.size();
    //int ncol=jmax-j0;
    // ==== lines
    glBegin(GL_LINES);
    glColor3f(0,0,0);
    int t=0;
    //printf(  "TableView Render %i %i %i %i \n", i0, j0, imax, jmax );
    t=ymin; for(int i=i0; i<imax;i++){ glVertex3f(xmin,t,0); glVertex3f(xmax,t,0); t+=fontSizeDef*2; }
    t=xmin; for(int j=j0; j<jmax;j++){ glVertex3f(t,ymin,0); glVertex3f(t,ymax,0); t+=fontSizeDef*nchs[j-j0]; }
    glEnd();
    //
    glColor3f(0,1.0,0);
    //int x=fontSizeDef*2;
    //int y=fontSizeDef*2;
    if( (i>=i0)&&(i<imax)&&(j>=j0)&&(j<jmax) ){
        Draw2D::drawRectangle( xmin+xs[j]-xs[j0], ymax+(i0-i)*2*fontSizeDef, xmin+xs[j+1]-xs[j0], ymax+(i0-i-1)*2*fontSizeDef, false );
    }
    // ==== text
    Draw  ::setRGB( textColor );
    char stmp[1024];
    for(int i=i0; i<imax;i++){
        int ch0 = 0;
        for(int j=j0; j<jmax;j++){
            int nch = table->toStr(i,j,stmp)-stmp;
            Draw2D::drawText( stmp, nch, {xmin+ch0*fontSizeDef, ymax-(i-i0+1)*fontSizeDef*2}, 0.0,  GUI_fontTex, fontSizeDef );
            ch0+=nchs[j];
        }
    }
}

void TableView::onKeyDown( const SDL_Event& e, GUI& gui ){
    if(input){
        input->onKeyDown( e );
    }
}

void TableView::onText( const SDL_Event& e, GUI& gui ){
    //printf( "TableView::onText() \n" );
    if( SDL_GetModState() & KMOD_CTRL ) return;
    if(input){
        input->inputText.insert(input->curPos,e.text.text); input->curPos++;
        printf( "TableView::onText() inputText >>%s<<\n", input->inputText.c_str() );
        input->modified = true;
    }
}

GUIAbstractPanel* TableView::onMouse( int x, int y, const SDL_Event& event, GUI& gui ){
    if( check( x, y ) ){
        //printf( "DropDownList::onMouse inside \n" );
        //if( event.type == SDL_MOUSEBUTTONUP ){
        if( event.type == SDL_MOUSEBUTTONDOWN ){
            int j_,i_ = i0 + (ymax-y)/(fontSizeDef*2);
            int dx=(x-xmin)+xs[j0];
            for(int j=j0;j<jmax;j++){ if(dx<xs[j+1]){ j_=j; break; } }
            if(event.button.button == SDL_BUTTON_LEFT){
                i=i_; j=j_;
                //printf( "TableView mouse select i,j %i %i\n", i, j );
                gui.bKeyEvents = false;
                SDL_StartTextInput();
                redraw = true;
            }
        }else if( event.type == SDL_MOUSEWHEEL ){
            //printf( " SDL_MOUSEWHEEL \n" );
            int ni = imax-i0;
            if     (event.wheel.y < 0){ i0 = _min( i0+1, table->n-ni ); }
            else if(event.wheel.y > 0){ i0 = _max( i0-1, 0           ); }
            imax=i0+ni;
            redraw = true;
        }
        return this;
    }
    return 0;
}

// ==============================
//    class GUI
// ==============================

GUIAbstractPanel* GUI::addPanel( GUIAbstractPanel* panel ){ panels.push_back(panel); return panels.back(); }

GUIAbstractPanel* GUI::onEvent( int mouseX, int mouseY, const SDL_Event& event ){
    GUIAbstractPanel* active = 0;
    //printf("GUI::onEvent \n");
    switch( event.type ){
        case SDL_KEYDOWN:
            //if(focused){ focused->onKeyDown( event ); }else{ txt.onKeyDown(  event ); }; break;
            if( focused && ( (event.key.keysym.sym == SDLK_TAB)||(event.key.keysym.sym == SDLK_RETURN)||(event.key.keysym.sym == SDLK_KP_ENTER) ) ){ 
                bTextEvents=!bTextEvents; 
                if(bTextEvents){ SDL_StartTextInput(); }else{ SDL_StopTextInput(); }
            }
            if(focused && bKeyEvents ){ 
                focused->onKeyDown( event, *this ); 
                active=focused; 
            }break;
        case SDL_TEXTINPUT:
            //if(focused){ focused->onText   ( event ); }else{ txt.onText   ( event );  }; break;
            if(focused && bTextEvents ){ 
                //printf( "GUI::onEvent() -> onText()  focused= `%s`| %li \n", focused->caption.c_str(), (long)focused );
                focused->onText( event, *this ); 
                active=focused;
            }break;
        case SDL_MOUSEWHEEL:
        case SDL_MOUSEBUTTONDOWN:
            active = NULL; focused=NULL;
            for(GUIAbstractPanel* panel: panels){
                active =  panel->onMouse( mouseX, mouseY, event, *this );
                if(active)focused=active;
            }break;
        case SDL_MOUSEBUTTONUP:
            if(event.button.button == SDL_BUTTON_LEFT){
                dragged = 0;
            }break;
        case SDL_MOUSEMOTION:
            SDL_MouseMotionEvent* event_ = (SDL_MouseMotionEvent*)&event;
            //if(GUI_mouse_panel) GUI_mouse_panel->moveTo( GUI_mouse_panel->xmin+event->xrel, GUI_mouse_panel->ymin+event->yrel );
            if(dragged){
                //printf(" GUI_globalEventHandler  SDL_MOUSEMOTION  %i %i \n", event_->xrel, -event_->yrel );
                dragged->moveBy( event_->xrel, -event_->yrel );
            }break;
    };
    return active;
}

void GUI::draw(){
    //printf( "GUI::draw() npanels=%i \n", panels.size() );
    //glLineWidth(1.0);
    //glLineWidth(0.5);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    for(GUIAbstractPanel* panel: panels){ if(focused!=panel)panel->draw(); }
    if(focused){
        if(focused->visible){ 
            focused->draw();
            if(bTextEvents){ glColor3f(1.0f,0.0f,0.0f); }else{Draw::setRGB(focused->textColor); }
            Draw2D::drawRectangle(focused->xmin,focused->ymin,focused->xmax,focused->ymax,false);
        }
    }
    //printf( "GUI::draw() END \n" );
}

void GUI::layoutRow( int xmin, int ymin, int xspace ){
    int x = xmin;
    for( GUIAbstractPanel* panel : panels ){
        int xsz = panel->xmax - panel->xmin;
        int ysz = panel->ymax - panel->ymin;
        panel->xmin = x;
        panel->xmax = x + xsz;
        panel->ymin = ymin;
        panel->ymax = ymin + ysz;
        x += xsz + xspace;
    };
}



