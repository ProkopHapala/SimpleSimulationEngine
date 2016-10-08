
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include <string>


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

	virtual void onKeyDown( SDL_Event e )                 = 0;
	virtual bool onMouse( int x, int y, SDL_Event event ) = 0;
    virtual void onText( SDL_Event e )                    = 0;

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
            Draw2D::setColor( bgColor );
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

	virtual void tryRender(){
		//Draw3D::drawRect( xmin, ymin, xmax, ymax );
		if(!redraw) return;
        gllist=glGenLists(1);
        glNewList( gllist, GL_COMPILE );
            glDisable( GL_LIGHTING );
            glDisable( GL_DEPTH_TEST);
            glShadeModel( GL_FLAT     );
            //glColor3f   ( 0.8, 0.8, 0.8 );
            Draw2D::setColor( bgColor );
            //printf("panel render %3.3f %3.3f %3.3f %3.3f \n", (float)xmin, (float)ymin, (float)xmax, (float)ymax );
            Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
            //float val_=(float)( (value-vmin)/(vmax-vmin) );
            if(isSlider){ Draw2D::setColor(barColor); Draw2D::drawRectangle ( xmin, ymax-2, xmin+val2x(value), ymax, true ); }
            Draw2D::setColor( bgColor );
            Draw2D::drawString ( caption,             xmin, ymin+12, 6, fontTex );
            //Draw2D::drawString ( val_text, 0, curPos, xmin, ymin,    6, fontTex );
            int nch = inputText.length();
            if( nch > 0 ){
                Draw2D::drawString ( inputText.c_str(), 0, nch, xmin, ymin,    6, fontTex );
                int xcur = xmin + curPos*6;
                Draw2D::setColor(textColor); Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+12} );
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
                    if(curPos>0) curPos--; redraw = true; break;
                case SDLK_RIGHT:
                    if(curPos<(inputText.length())) curPos++; redraw = true; break;
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
	}

	virtual bool onMouse( int x, int y, SDL_Event event ){
        //printf( "panel.onMouse %i %i \n", x, y );
		if( check( x, y ) ){
			toRelative(x,y);
			//printf( "  panel.onMouse %i %i \n", x, y );
			if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
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
		return executed;
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
        glPushMatrix();
        glTranslatef(xmin,ymin,0);
        for(int i=0; i<nsubs; i++){
            subs[i]->draw();
        }
        glPopMatrix();
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
            subs[i]->init(0,yi,xmax-xmin,yi+dy, fontTex );
            subs[i]->caption = new char[16];
            sprintf(subs[i]->caption,"val%i",i);
            yi+=dy;
        }
	};

	virtual bool onMouse  ( int x, int y, SDL_Event event ){
        if( check( x, y ) ){
			toRelative(x,y);
            for(int i=0; i<nsubs; i++){
                subs[i]->onMouse  ( x, y, event );
            }
        }
	};

    virtual void onKeyDown( SDL_Event e ){};

    virtual void onText   ( SDL_Event e ){};

};




#endif
