
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include "Table.h"

#include <string>

/*

TODO:
 - draw small glyph by simple language ... like:
    goto    100.0 100.0
    line    200.0 250.0
    ngon    6 16.5
*/

//class GUIAbstractPanel;

//extern GUIAbstractPanel* GUI_mouse_panel;
//extern Vec2i GUI_mouse_old_pos;
extern int GUI_fontTex;

//void GUI_globalEventHandler(const SDL_Event* event );

class GUI;


class GUIAbstractPanel;

// Alterntively we may use events
//    https://wiki.libsdl.org/SDL_UserEvent
//    https://stackoverflow.com/questions/26707536/event-driven-programming-callback-vs-message-polling
class GUIEventCallback{public:
    virtual int GUIcallback(GUIAbstractPanel* caller)=0;
};

// ==============================
//    class GUITextInput
// ==============================

class GUITextInput{ public:

    int         curPos=0;
	std::string inputText;

    bool     isNumber=true,checkRange=false;
	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;
	SDL_Keycode num_op = 0;

	bool     modified=true,entered=false;

	// ==== functions

	void         applyVal( float f );

    virtual void view3D( const Vec3d& pos, int fontTex, float textSize );
    virtual void viewHUD( const Vec2i& pos, int fontTex );
    virtual void onKeyDown( SDL_Event e );
	virtual void onText   ( SDL_Event e );

};

// ==============================
//    class GUIAbstractPanel
// ==============================

class GUIAbstractPanel{ public:
    //int textSz = fontSizeDef;
	int  xmin=256,xmax=128,ymin=0,ymax=0;
	bool visible=true, disabled=false;

	uint32_t bgColor=0xA0A0A0, textColor=0x000000;

	bool     redraw=true;
	int      gllist=0;

	//int      fontTex=0;
    //char*    caption=NULL;
    std::string caption;

	// ==== functions

    void initPanel( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ );

    GUIAbstractPanel(){};
    GUIAbstractPanel( const std::string& caption, int xmin, int ymin, int xmax, int ymax ){ initPanel(caption, xmin,ymin,xmax,ymax); };
    virtual ~GUIAbstractPanel(){};

	virtual void moveTo(int x, int y);
	virtual void moveBy(int dx,int dy);

	virtual void              onKeyDown( const SDL_Event& e, GUI& gui );
	virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui );
    virtual void              onText( const SDL_Event& e );

    virtual void view  ( );
    //virtual void tryRender();
    //void view     ();
    void tryRender();
    virtual void render();

    // inline fnctions

    inline void draw      ( ){ tryRender(); view(); };
    inline bool check     ( int  x, int  y ){
        //printf( "check x %i <%i...%i>   y %i <%i...%i>\n", x, xmin, xmax,   y, ymin, ymax );
        return (x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax);
    }
	inline void toRelative ( int& x, int& y ){ x-=xmin; y-=ymin; }

};

// ==============================
//       class GUIPanel
// ==============================

class GUIPanel : public GUIAbstractPanel { public:
	bool isSlider=true, isButton=false;

	uint32_t barColor=0x00FF00;

	bool     executed=false;
	int      curPos=0;
	std::string inputText;

	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0d;
	bool     isInt = false;

	void (*command)(double) = NULL;

    // ==== functions

    GUIPanel()=default;
    GUIPanel( const std::string& caption, int xmin, int ymin, int xmax, int ymax, bool isSlider_, bool isButton_ ){
        initPanel(caption, xmin,ymin,xmax,ymax); isSlider=isSlider_; isButton=isButton_;
        command=0;
    };

    virtual void view();
	//virtual void tryRender();
	void render();
    virtual void              onKeyDown( const SDL_Event&  e );
    virtual void              onText( const SDL_Event&  e );
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui );

	// ===== inline functions
    inline  void  val2text()         { inputText = std::to_string(value); };
	inline double x2val( float  x   ){ return ( x*(vmax-vmin)/(xmax-xmin) )+ vmin; };
	inline float  val2x( double val ){ return (val-vmin)*(xmax-xmin)/(vmax-vmin);  };

	inline GUIPanel* setRange(float vmin_, float vmax_){ vmin=vmin_; vmax=vmax_; return this; };

};

// ==============================
//     class  ScisorBox
// ==============================

class ScisorBox : public GUIAbstractPanel { public:

    // ==== functions

    void apply();
    void initScisor( const std::string& caption, int xmin_, int ymin_, int xmax_, int ymax_ );

    ScisorBox(){};
    ScisorBox( const std::string& caption, int xmin_, int ymin_, int xmax_, int ymax_ ){ initScisor( caption, xmin_, ymin_, xmax_, ymax_ ); };

    //virtual void draw     ( );
    //virtual void tryRender( );
    virtual void render( );
    virtual GUIAbstractPanel* onMouse ( int x, int y, const SDL_Event&  event, GUI& gui );

    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};

};

// ==============================
//     class  MultiPanel
// ==============================

class MultiPanel : public GUIAbstractPanel { public:
    int nsubs;
    GUIPanel ** subs;
    bool opened = true;
    int dy;

    // ==== functions

    void initMulti( const std::string& caption, int xmin_, int ymin_, int xmax_, int dy, int nsubs_ );

    MultiPanel(){};
    MultiPanel(const std::string& caption, int xmin, int ymin, int xmax, int dy, int nsubs){ initMulti( caption, xmin, ymin, xmax, dy, nsubs ); }

    virtual void open();
    virtual void close();
    void toggleOpen();

    virtual void moveBy(int dx, int dy);

    virtual void view  ( );
    //virtual void tryRender( );
    virtual void render( );
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui );

    //virtual void onKeyDown( SDL_Event e, GUI& gui ){};
    //virtual void onText   ( SDL_Event e, GUI& gui ){};

};

// ==============================
//     class  DropDownList
// ==============================

//static const char* exampleDropDownListItems[3] = {"Item1","Item2","Item3"};

class DropDownList : public GUIAbstractPanel { public:
    bool bOpened = true;
    int nSlots=5;
    //int iSlot0=0;

    int iSelected=0;
    int iItem0 = 0;
    //int nItems = 0;
    //char** labels = (char**)exampleDropDownListItems;
    std::vector<std::string> labels;

    GUIEventCallback* onSelect = 0;

    //int nsubs;
    //GUIPanel ** subs;

    // ==== functions

    DropDownList* addItem(const std::string& label);
    void initList( const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_ );
    int selectedToStr(char* str);

    DropDownList(){}
    DropDownList( const std::string& caption, int xmin, int ymin, int xmax, int nSlots){ initList(caption,xmin,ymin,xmax,nSlots); }

    virtual void open();
    virtual void close();

    //virtual void view ( );
    //virtual void tryRender( );
    virtual void render( );
    virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui );

    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};

};



// ==============================
//     class  TreeView
// ==============================

#include "Tree.h"

class TreeViewItem{ public:
    //bool open=false;
    bool open=true;
    int level=0;
    int nth_line=0;
    std::string caption;

    TreeViewItem(){};
    TreeViewItem( std::string caption_){ caption=caption_; };
};
typedef PTree<TreeViewItem> TreeViewTree;

//static const char* exampleDropDownListItems[3] = {"Item1","Item2","Item3"};
class TreeView : public GUIAbstractPanel { public:
    //typedef Tree<TreeViewItem> TreeT;
    int iItem0=0;
    int iSelected=0;
    int nSlots = 5;
    TreeViewTree root;
    std::vector<TreeViewTree*> lines;


    virtual void view();
    virtual void render();
    void updateLines( TreeViewTree& node, int level );

    inline void updateLines(){ lines.clear(); updateLines(root,0); };

    void initTreeView( const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_ );

    TreeView(){}
    TreeView( const std::string& caption, int xmin, int ymin, int xmax, int nSlots){
        initTreeView(caption,xmin,ymin,xmax,nSlots);
    }

    virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui );
    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};
};

// ==============================
//     TableView
// ==============================

//static const char* exampleDropDownListItems[3] = {"Item1","Item2","Item3"};
class TableView : public GUIAbstractPanel { public:
    //typedef Tree<TreeViewItem> TreeT;
    Table* table=0;
    int i0=0,j0=0;
    int imax=0,jmax=0;
    int i=0,j=0;  // cursors
    std::vector<int> nchs;
    std::vector<int> xs;

    //virtual void view();
    //virtual void render();
    //void updateLines( TreeViewTree& node, int level );
    //inline void updateLines(){ lines.clear(); updateLines(root,0); };

    void initTableView( Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_ ){
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

        printf( " i (%i,%i) j (%i,%i) \n", i0, imax,    j0, jmax   );
        printf( " x (%i,%i) y (%i,%i) \n", xmin, xmax,  ymin, ymax );

        redraw = true;
    };

    TableView()=default;
    TableView( Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_ ){
        initTableView( table_, caption_, xmin_, ymin_, i0_, j0_, imax_, jmax_ );
    }
    //TreeView( const std::string& caption, int xmin, int ymin, int xmax, int nSlots){
    //    initTreeView(caption,xmin,ymin,xmax,nSlots);
    //}

    //virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui );
    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};

    // Implementaation

    //virtual void view();
    virtual void render(){
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
        printf(  "TableView Render %i %i %i %i \n", i0, j0, imax, jmax );
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

    void view ( ){
        // Draw2D::drawPointCross({xmin,ymax},5);
        glCallList( gllist );
        //int xcur = xmin + curPos*fontSizeDef;
        //Draw2D::drawLine   ( {xcur, ymin}, {xcur, ymin+fontSizeDef*2} );
    };

    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui ){
        if( check( x, y ) ){
            //printf( "DropDownList::onMouse inside \n" );
            //if( event.type == SDL_MOUSEBUTTONUP ){
            if( event.type == SDL_MOUSEBUTTONDOWN ){
                int j_,i_ = i0 + (ymax-y)/(fontSizeDef*2);
                int dx=(x-xmin)+xs[j0];
                for(int j=j0;j<jmax;j++){ if(dx<xs[j+1]){ j_=j; break; } }
                if(event.button.button == SDL_BUTTON_LEFT){
                    i=i_; j=j_;
                    printf( "TableView mouse select i,j %i %i\n", i, j );
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
    };

};


// ==============================
//    class GUI
// ==============================

class GUI{ public:
    GUIAbstractPanel* focused = 0;
    GUIAbstractPanel* dragged = 0;
    std::vector<GUIAbstractPanel*> panels;

    GUIAbstractPanel* addPanel( GUIAbstractPanel* panel );
    GUIAbstractPanel* onEvent( int mouseX, int mouseY, const SDL_Event& event );
    void draw();

    ~GUI(){
        for(GUIAbstractPanel* panel: panels){ delete panel; }
    }

};

#endif
