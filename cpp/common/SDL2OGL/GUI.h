
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

class GUIPanel{
    public:
	// ===== properties
	int  xmin=256,xmax=128,ymin=0,ymax=0;
	bool visible=true, disabled=false;
	bool isSlider=true, isButton=false;

	uint32_t bgColor=0x000000, textColor=0xFFFFFF, barColor=0x00FF00;
	int      fontTex=0;
	bool     redraw=true;
	bool     executed=false;
	int      gllist=0;

    //static const int nSubMax=64;
	//int       nSub;
	//GUIPanel  subs[nSubMax];

	int      curPos=0;
	int      nCharMax=64;
	char*    val_text=NULL;
	char*    caption =NULL;

	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;

	void (*command)(double) = NULL;

	// ===== inline functions
	inline bool check      ( int  x, int  y ){  return (x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax); }
	inline void toRelative ( int& x, int& y ){ x-=xmin; y-=ymin; }

	// ===== virtual functions

	virtual void render(){
		//Draw3D::drawRect( xmin, ymin, xmax, ymax );
       // gllist=glGenLists(1);
       // glNewList( gllist, GL_COMPILE );
            glDisable( GL_LIGHTING );
            glDisable( GL_DEPTH_TEST);
            glShadeModel( GL_FLAT     );
            //glColor3f   ( 0.8, 0.8, 0.8 );
            Draw2D::setColor( bgColor );
            //printf("panel render %3.3f %3.3f %3.3f %3.3f \n", (float)xmin, (float)ymin, (float)xmax, (float)ymax );
            Draw2D::drawRectangle ( xmin, ymin, xmax, ymax, true );
            float val_=(float)( (value-vmin)/(vmax-vmin) );
            if(isSlider){ Draw2D::setColor(barColor); Draw2D::drawRectangle ( xmin, ymin, xmin+val_*(xmax-xmin), ymax, true ); }
            Draw2D::drawString ( caption,  xmin, ymin+12, 6, fontTex );
            Draw2D::drawString ( val_text, xmin, ymin,    6, fontTex );
       // glEndList();
	};

	/*
	virtual void onEvent( SDL_Event e ){
		// see https://wiki.libsdl.org/SDL_Keysym
		if( e.type == SDL_KEYDOWN ){
			if( e.key.keysym.sym == SDLK_BACKSPACE && inputText.length() > 0 )	{
				inputText.pop_back();
				renderText = true;
			}else if( e.key.keysym.sym == SDLK_c && SDL_GetModState() & KMOD_CTRL )	{
				SDL_SetClipboardText( inputText.c_str() );
			}else if( e.key.keysym.sym == SDLK_v && SDL_GetModState() & KMOD_CTRL )	{
				inputText = SDL_GetClipboardText();
				renderText = true;
			}else if( (e.key.keysym.sym == SDLK_RETURN) || ( e.key.keysym.sym == SDLK_KP_ENTER) ){
				printf( " lua do : %s \n", inputText.c_str() );
				luaL_dostring(L, inputText.c_str() );
				inputText = " ";
				renderText = true;
			}
		}else if( e.type == SDL_TEXTINPUT )	{
			//if( !( ( e.text.text[ 0 ] == 'c' || e.text.text[ 0 ] == 'C' ) && ( e.text.text[ 0 ] == 'v' || e.text.text[ 0 ] == 'V' ) && SDL_GetModState() & KMOD_CTRL ) ){
				inputText += e.text.text;
				try{ val_float = std::stof( inputText ); }catch (int e){ isFloat=false; isInt=false; }
				redraw     = true;
			//}
		}
	}
   */

	virtual bool onMouse( int x, int y, SDL_Event event ){
        //printf( "panel.onMouse %i %i \n", x, y );
		if( check( x, y ) ){
			toRelative(x,y);
			//printf( "  panel.onMouse %i %i \n", x, y );
			if( ( event.type == SDL_MOUSEBUTTONDOWN ) ){
                if(isSlider && (event.button.button==SDL_BUTTON_RIGHT)){
                    value=( x*(vmax-vmin)/(xmax-xmin) ) + vmin;
                    sprintf(val_text, "%3.3f", value );
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

	void init( int xmin_, int ymin_, int xmax_, int ymax_, int fontTex_ ){
		xmin=xmin_,ymin=ymin_,xmax=xmax_,ymax=ymax_; fontTex=fontTex_;
		val_text = new char[nCharMax];
	}

};

#endif
