
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include <string>

// ==============================
//    class GUIPanelBasic
// ==============================

class GUITextInput{
    public:

    int         curPos=0;
	std::string inputText;

    bool     isNumber=true,checkRange=false;
	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;

	bool     modified=true,entered=false;

    SDL_Keycode num_op = 0;
    //float    textSize  = 6.0;
	//uint32_t textColor = 0x000000;

	void applyVal( float f ){
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

    virtual void view3D( const Vec3d& pos, int fontTex, float textSize ){
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

    virtual void onKeyDown( SDL_Event e ){
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


	virtual void onText( SDL_Event e ){
        if( SDL_GetModState() & KMOD_CTRL ) return;
        inputText.insert(curPos,e.text.text); curPos++;
        modified = true;
	}

};

// ==============================
//    class GUIPanelBasic
// ==============================

class GUIAbstractPanel{
    public:
	int  xmin=256,xmax=128,ymin=0,ymax=0;
	bool visible=true, disabled=false;

	uint32_t bgColor=0xA0A0A0, textColor=0x000000;

	bool     redraw=true;
	int      gllist=0;

	int      fontTex=0;
    char*    caption=NULL;

    inline bool check      ( int  x, int  y ){  return (x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax); }
	inline void toRelative ( int& x, int& y ){ x-=xmin; y-=ymin; }

	virtual void              onKeyDown( SDL_Event e )                 = 0;
	virtual GUIAbstractPanel* onMouse( int x, int y, SDL_Event event ) = 0;
    virtual void              onText( SDL_Event e )                    = 0;

    virtual void draw  ( ){
        glCallList( gllist );
    };

    virtual void tryRender(){
        if(!redraw) return;
        gllist=glGenLists(1);
        glNewList( gllist, GL_COMPILE );
            glDisable   ( GL_LIGHTING    );
            glDisable   ( GL_DEPTH_TEST  );
            glShadeModel( GL_FLAT        );
            Draw  ::setRGB( bgColor );
            Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
            if(caption) Draw2D::drawString ( caption,             xmin, ymax-12, 6, fontTex );
        glEndList();
        redraw=false;
	};

    virtual void init( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_ ){
		xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
		//val_text = new char[nCharMax];
	};

};


// ==============================
//       GUIPanel
// ==============================


class GUIPanel : public GUIAbstractPanel {
    public:
	// ===== properties

	bool isSlider=true, isButton=false;

	uint32_t barColor=0x00FF00;

	bool     executed=false;

    //static const int nSubMax=64;
	//int       nSub;
	//GUIPanel  subs[nSubMax];

	int      curPos=0;
	//int      nCharMax=64;
	//char*    val_text=NULL;
	std::string inputText;

	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;

	void (*command)(double) = NULL;

	// ===== inline functions
	inline double x2val( float  x   ){ return ( x*(vmax-vmin)/(xmax-xmin) )+ vmin; };
	inline float  val2x( double val ){ return (val-vmin)*(xmax-xmin)/(vmax-vmin);  };

	// ===== virtual functions

    virtual void draw  ( ){
        glCallList( gllist );
        int xcur = xmin + curPos*6;
        Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+12} );
    };


	virtual void tryRender(){
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

	virtual void onKeyDown( SDL_Event e ){
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


	virtual void onText( SDL_Event e ){
        if( SDL_GetModState() & KMOD_CTRL ) return;
        //char ch = e.text.text[0];
        //printf( "input event >>%s<<\n", e.text.text );
        //if((curPos+1)<nCharMax)curPos++;
        //val_text[curPos] = ch;
        //inputText.push_back(ch);
        inputText.insert(curPos,e.text.text); curPos++;
        redraw = true;
	}

	virtual GUIAbstractPanel* onMouse( int x, int y, SDL_Event event ){
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

};

// ==============================
//       GUIPanel
// ==============================

class MultiPanel : public GUIAbstractPanel {
    public:
    int nsubs;
    GUIPanel ** subs;

    virtual void draw  ( ){
        glCallList( gllist );
        for(int i=0; i<nsubs; i++){
            subs[i]->draw();
        }
	};

    virtual void tryRender( ){
        GUIAbstractPanel::tryRender();
        for(int i=0; i<nsubs; i++){
            subs[i]->tryRender();
        }
	};

    void initMulti( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_, int nsubs_ ){
		xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
		nsubs=nsubs_;
		subs = new GUIPanel*[nsubs_];
		int dy = (ymax-ymin)/nsubs;
		int yi = 0;
        for(int i=0; i<nsubs; i++){
            subs[i] = new GUIPanel();
            subs[i]->init(xmin,yi,xmax,yi+dy, fontTex );
            subs[i]->caption = new char[16];
            sprintf(subs[i]->caption,"val%i",i);
            yi+=dy;
        }
	};

	virtual GUIAbstractPanel* onMouse  ( int x, int y, SDL_Event event ){
        GUIAbstractPanel* active = NULL;
        if( check( x, y ) ){
            for(int i=0; i<nsubs; i++){
                active = subs[i]->onMouse ( x, y, event );
                if(active) return active;
            }
        }
        return active;
	};

    virtual void onKeyDown( SDL_Event e ){};
    virtual void onText   ( SDL_Event e ){};

};




#endif
