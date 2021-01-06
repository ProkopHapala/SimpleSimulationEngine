
#ifndef  ShapePainter_h
#define  ShapePainter_h

#include "Draw.h"
#include "Draw2D.h"

class Brush{ public:
    virtual void paint()=0;
    virtual void view ()=0;
    virtual void eventHandling( const SDL_Event& event, Vec2f mpos, uint32_t mouseState )=0;
};


// ==========================
//      CornerBrush
// ==========================

class CornerBrush : public Brush{ public:

uint32_t color = 0xFF808080;
Vec2f pos;   // position of spike
Vec2f pend;  // position of 1st base point
Vec2f rot;   // rotation of second end with respect to first

inline void getRot( const Vec2f& p ){
    rot.set_div_cmplx( p, pend-pos );
    rot.normalize();
}

inline Vec2f getC(){ Vec2f C; C.set_mul_cmplx(pend-pos, rot); return pos+C; }

/*
void normalize(){
    Vec2f d1 = pend-pos;
    float r2 = d1.norm2();
    d1 = pC-pos;
    d1.mul( sqrt(r2/d1.norm2()) );
    pC = pos + d1;
}
*/

virtual void paint()override{  // deposit on canvas
    Vec2f pC = getC();
    glBegin(GL_TRIANGLES);
    Draw::setRGBA( color             ); glVertex2f( pos.x , pos.y  );
    Draw::setRGBA( color &0x00FFFFFF ); glVertex2f( pend.x, pend.y );
    Draw::setRGBA( color &0x00FFFFFF ); glVertex2f( pC.x  , pC.y   );
    glEnd();
}

virtual void view()override{ // show on canvas without depositing
    //printf( " (%g,%g)  (%g,%g) (%g,%g) \n", pos.x,pos.y,  pend.x,pend.y, pC.x,pC.y );
    Vec2f pC = getC();
    float psz=0.1;
    glColor3f( 1.0,0.0,0.0 ); Draw2D::drawPointCross( pos,  psz );
    glColor3f( 0.0,1.0,0.0 ); Draw2D::drawPointCross( pend, psz );
    glColor3f( 0.0,0.0,1.0 ); Draw2D::drawPointCross( pC,   psz );
    glBegin(GL_TRIANGLES);
    //glBegin(GL_LINE_LOOP);
    Draw::setRGBA( color             ); glVertex2f( pos.x , pos.y  );
    Draw::setRGBA( color &0x00FFFFFF ); glVertex2f( pend.x, pend.y );
    Draw::setRGBA( color &0x00FFFFFF ); glVertex2f( pC.x  , pC.y   );
    glEnd();
}

virtual void eventHandling ( const SDL_Event& event, Vec2f mpos, uint32_t mouseState )override{
    //printf(  "mouseState : LMB %i  RMB %i \n", mouseState & SDL_BUTTON(SDL_BUTTON_LEFT),   mouseState & SDL_BUTTON(SDL_BUTTON_RIGHT)    );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_m:  edit_mode = (EDIT_MODE)((((int)edit_mode)+1)%((int)EDIT_MODE::size)); printf("edit_mode %i\n", (int)edit_mode); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_l:    onSelectLuaShipScript.GUIcallback(lstLuaFiles); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  pos=mpos;  break;
                //case SDL_BUTTON_RIGHT: pC =mpos;  break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT: pend=mpos; paint(); break;
                //case SDL_BUTTON_RIGHT:break;
            }
            break;
        case SDL_MOUSEMOTION:
            //case SDL_BUTTON_LEFT:  pend=mpos; normalize(); break;
            //case SDL_BUTTON_RIGHT: pC  =mpos; normalize(); break;
            if(  mouseState & SDL_BUTTON(SDL_BUTTON_LEFT ) ){ pend=mpos;       }
            if(  mouseState & SDL_BUTTON(SDL_BUTTON_RIGHT) ){ getRot( mpos );  }
        break;
    };
}

};

#endif
