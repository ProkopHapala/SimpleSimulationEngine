
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>



#include "drawMath2D.h"
//#include "drawDrawUtils.h"
#include "Voronoi.h"

#include "GameScreen.h" // THE HEADER





void GameScreen::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( " GameScreen::draw \n" );

	glEnable (GL_LIGHTING);
	glShadeModel(GL_FLAT);

	world->update();

    //printf( " world->update DONE \n" );

	//for( std::vector<Projectile*>::iterator it = world->projectiles.begin(); it != world->projectiles.end(); ++it ) {
	for( auto ship : world->ships ) {
        //printf( "projectile draw \n " );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->drawHitBox( );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->draw_shape( );
        glColor3f( 0.2f, 0.2f, 0.2f );  ship->draw( );
	}

	for( auto p : world->projectiles ) {
        //printf( "projectile draw \n " );
		p->draw();
		//p->update_old_pos();
	}

	// ============== Convex Polygons Mouse Raycast

/*
    //printf( " mouse x,y : %i %i %f %f \n", mouseX, mouseY, mouseRight(mouseX), mouseUp(mouseY) );
    Vec2d mvec;
    mvec.set( mouseRight(mouseX), mouseUp(mouseY) );
    //mvec.set( 0.0, 0.0 );
    Draw2D::drawPointCross_d( mvec, 1.0 );

	for( auto isle : world->isles ) {
        //printf( "draw isle: %i %f %f \n ", isle->n, isle->corners[0].x, isle->corners[0].y );
        if( isle->pointIn( mvec ) ){
            glColor3f( 1.0f, 0.2f, 0.2f );
            //printf( " mouse in isle \n " );
        }else{
            glColor3f( 0.2f, 0.2f, 0.2f );
        }
		Draw2D::drawConvexPolygon( isle->n, isle->corners, true );

		//exit(0);
	}
*/

	// ============== Convex Polygons Cut

/*

    Convex2d base( 5 );

	//base.corners[0].set( -1.0, +1.0 );
    //base.corners[1].set( +1.0, -1.5 );
    //base.corners[2].set( +2.0,  0.5 );

	base.corners[0].set( -1.0, -1.0 );
    base.corners[1].set( +1.0, -1.0 );
    base.corners[2].set( +2.0,  0.0 );
	base.corners[3].set( +1.0, +1.0 );
	base.corners[4].set( -1.0, +1.0 );

    glColor3f( 0.9f, 0.9f, 0.9f ); Draw2D::drawConvexPolygon( base.n, base.corners, false );


    Vec2d Acut,Bcut;
    Acut.set( -1.0, -5.0 );
    Bcut.set( +3.0, +5.0 );
    Line2d cutline; cutline.set( Acut, Bcut );
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine_d( Acut, Bcut );

    Convex2d left;
    Convex2d right;

    base.cut( cutline, left, right );
    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::drawPointCross_d( left.corners[0], 0.3 );
    glColor3f( 0.2f, 0.9f, 0.9f ); Draw2D::drawPointCross_d( left.corners[1], 0.3 );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawConvexPolygon( left.n,  left.corners,  true );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::drawConvexPolygon( right.n, right.corners, true );

*/

    // ============== plot Voronoi

/*
    glColor3f( 0.9f, 0.9f, 0.9f );
    int iii=0;
    for ( Vec2d * p : *(world->voronoi->places) ){
        iii++;
        Draw2D::drawPointCross_d( *p, 0.1 );
	}
	//printf( " Vertices N %i \n", iii );

    glColor3f( 0.9f, 0.9f, 0.9f );
    iii=0;
    for ( VoronoiNamespace::VEdge * pe : *(world->voronoi->edges) ){
        iii++;
        Vec2d a,b;
        a.set( *(pe->start) );
        b.set( *(pe->end)   );
        Draw2D::drawLine_d( a, b );
        //printf( " %f %f   %f %f \n", a.x, a.y,    b.x, b.y  );
	}
	//exit(0);
	//printf( " Edges N %i \n", iii );

*/

    // ============== polar

    if(thisShip!=NULL){
/*
        const int nsteps = 100;
        Vec2d * vels = new Vec2d[ nsteps ];
        Vec2d * rots = new Vec2d[ nsteps ];
        Vec2d * poss = new Vec2d[ nsteps ];
        thisShip->pos.set( 0.0d );
        thisShip->vel.set( 0.0d );
        glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );
        thisShip->testSail( nsteps, 20, 0.01, -10.0, poss, vels, rots );
        Draw2D::drawLines( nsteps, poss );
        delete vels,rots,poss;
*/

/*
        Vec2d vel_conv, rot_conv;
        int nstepmax  = 1000;
        thisShip->pos.set( 0.0d );
        thisShip->vel.set( 0.0d );
        glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );
        int nstepconv = thisShip->convergeSail( nstepmax, 20, 0.001, -10.0, 1e-3, 1e-3, vel_conv, rot_conv );
        glColor3f( 0.9f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );
        if ( nstepconv < nstepmax ){
            printf( " convergeSail : %i (%3.3f,%3.3f) (%3.3f,%3.3f) \n", nstepconv, vel_conv, rot_conv );
        }else{
            printf( " convergeSail not converged int %i ! \n", nstepmax );
        }
*/

        const int nsteps = 50;
        double * phi_rudder = new double[ nsteps ];
        double * phi_mast   = new double[ nsteps ];
        double * wind_speed = new double[ nsteps ];
        Vec2d  * vels        = new Vec2d[ nsteps ];
        Vec2d  * rots        = new Vec2d[ nsteps ];
        glColor3f( 0.2f, 0.2f, 0.2f );
        for( int i=0; i<nsteps; i++ ){
            //phi_rudder[ i ] = 1.57079632679 +  -0.3 + ( i*( 0.3-(-0.3) )/ nsteps );
            //phi_mast  [ i ] = M_PI * 0.25;
            //phi_mast  [ i ] = 0.0;

            phi_rudder[ i ] = 1.57079632679 + 0.09;
            phi_mast  [ i ] = M_PI * (  -0.5   +    i/(float)nsteps  );
            phi_mast  [ i ] = M_PI * ( 0.0   +  i/(float)nsteps  );


            //phi_rudder[ i ] = 1.57079632679 + 0.1 * i/(float)nsteps ;
            //phi_mast  [ i ] = M_PI * +0.25;

            wind_speed[ i ] = -10.0d;
            vels[ i ].set( 0.0d, 0.0d );
            rots[ i ].set( 1.0d, 0.0d ); rots[ i ].normalize();
        }
        thisShip->evalPolar( nsteps, 0.01, 1e-300, 1e-300,  phi_rudder, phi_mast, wind_speed, vels, rots, true );
        glScalef( 10.0, 10.0, 10.0 );
        glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawLines( nsteps, vels );
        glScalef( 1/10.0, 1/10.0, 1/10.0 );
        delete phi_rudder, phi_mast, wind_speed,  vels, rots;


        STOP = true;

    }


	// ============== rest

	glDisable  (GL_LIGHTING);

	//exit(0);
	//drawAxis( 10 );
};


void GameScreen::drawHUD(){
    float bar_y =  10.0f;
    float bar_x = 100.0f;
    if( thisShip != NULL ){
        Draw2D::drawRectangle( { 0.0f,  0.0f }, { bar_x * thisShip->gunload_left,    bar_y }, true );
        Draw2D::drawRectangle( { 0.0f, bar_y }, { bar_x * thisShip->gunload_right, 2*bar_y }, true );
        Draw2D::drawLine     ( { bar_x, 0.0f }, { bar_x, 2*bar_y } );
    }

    Vec2d compass_pos; compass_pos.set( WIDTH*0.9, HEIGHT*0.9 );
	glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( compass_pos, 10 );
	glColor3f( 0.2f, 0.5f, 0.2f );  Draw2D::drawVecInPos_d  ( ( *(Vec2d*)&(world->wind_speed))*10.1,   compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  Draw2D::drawVecInPos_d  ( world->watter_speed*10.1, compass_pos );
};




GameScreen::GameScreen( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) {
//	init( id, WIDTH_, HEIGHT_ );
}

