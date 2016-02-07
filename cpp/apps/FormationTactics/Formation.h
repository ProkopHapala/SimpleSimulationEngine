

#ifndef Formation_h
#define Formation_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "drawMath2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Body2D.h"

#include "Soldier.h"

class Formation{
	public:
    char  * name;

    Vec2d * nodes;

    Vec2d center;
    Vec2d p00target,p01target;
    bool movingToTarget = false;
    //Vec2d pL,pR;
    Vec2d p00,p01,p10,p11;
    Vec2d dirFw,dirLf;
    double length,width;
    double kLength=3.0, kWidth=3.0;

    int nrows;
    int ncols;
    int nsoldiers;
    int ncapable,nalive;
    Soldier * soldiers = NULL;


    void moveBy( const Vec2d& dpos ){
        center.add( dpos );
        setCenterRot( center, dirFw );
    }


    void interactInside( ){
        for( int i=0; i<nsoldiers; i++ ){
            Soldier * si = soldiers + i;
            if( si->alive ){
                for( int j=0; j<i; j++ ){
                    Soldier * sj = soldiers + j;
                    if( sj->alive ){
                        Vec2d d;
                        d.set_sub( si->pos, sj->pos );
                        double r2 = d.norm2( );
                        if( r2 < 1.0 ){
                            d.mul( (1-r2)*10 );
                            si->force.add( d );
                            sj->force.sub( d );
                        }
                    }
                }
            }
        }
    }

    void setTarget( const Vec2d& target ){
        printf( " setTarget (%3.3f,%3.3f) \n", target.x, target.y );
        Vec2d d;
        d.set_sub( target, center ); d.normalize();
        //dirFw.set( d ); dirFw.normalize(); dirLf.set_perp( dirFw );
        setCenterRot( center, d );
        p00target.set_add_mul( target, dirLf,  length * 0.5d );
        p01target.set_add_mul( target, dirLf, -length * 0.5d );
        movingToTarget = true;
    }

    void checkTarget( ){
        Vec2d d;
        d.set_sub( p00, p00target );
        if( d.norm2() >  0.01 ) return;
        d.set_sub( p01, p01target );
        if( d.norm2() >  0.01 ) return;
        movingToTarget = false;
    }

    void moveToTarget( ){
        checkTarget( );
        if( movingToTarget ){
            double speed = 0.01;
            Vec2d d1,d0;
            d0.set_sub( p00target, p00 ); d0.normalize(); p00.add_mul( d0, speed );
            d1.set_sub( p01target, p01 ); d1.normalize(); p01.add_mul( d1, speed );
            setEnds( p00, p01, width );
        }
    }



    void applyWillForce( Soldier& soldier ){
        if( soldier.capable ){
            Vec2d d;
            d.set_sub( soldier.pos, p00 );
            double llf = -dirLf.dot( d );
            if      ( llf < 0      ){ soldier.willForce.add_mul( dirLf,  llf        *kLength ); }
            else if ( llf > length ){ soldier.willForce.add_mul( dirLf, (llf-length)*kLength ); }
            double lfw = -dirFw.dot( d );
            if      ( lfw < 0      ){ soldier.willForce.add_mul( dirFw,  lfw        *kWidth  ); }
            else if ( lfw > width  ){ soldier.willForce.add_mul( dirFw, (lfw-width )*kWidth  ); }
            //printf( "- %3.3f %3.3f   %3.3f %3.3f \n", llf, lfw, length, width  );
        }
    }

    void applyWillForce( ){
        for( int i = 0; i<nsoldiers; i++ ){
            applyWillForce( soldiers[i] );
            //printf( " (%3.3f,%3.3f)\n", i, soldiers[i].willForce.x, soldiers[i].willForce.y );
        }
    }

    void clean_temp(){
        for( int i = 0; i<nsoldiers; i++ ){
            soldiers[i].clean_temp();
        }
    }

    void setEnds( const Vec2d& pL_, const Vec2d& pR_, double width_ ){
        width = width_;
        p00.set( pL_ ); p01.set( pR_ );
        center.set_lincomb( 0.5d, p00, 0.5d, p01 );
        Vec2d d; d.set_sub( p00, p01 );
        length = d.norm();
        dirLf.set_mul ( d, 1/length  );
        dirFw.set_perp( dirLf );
        p10.set_add_mul( p00, dirFw, -width );
        p11.set_add_mul( p01, dirFw, -width );
        //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
    }


    void setCenterRot( const Vec2d& center_, const Vec2d& dirFw_ ){
        center.set( center_ );
        dirFw .set( dirFw_ );
        dirLf .set_perp ( dirFw );
        p00.set_add_mul( center, dirLf,  length * 0.5d );
        p01.set_add_mul( center, dirLf, -length * 0.5d );
        p10.set_add_mul( p00,    dirFw, -width );
        p11.set_add_mul( p01,    dirFw, -width );
        //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
        //exit(0);
    }

    void update( double dt ){
        moveToTarget( );
        for( int i = 0; i<nsoldiers; i++ ){
            soldiers[i].vel.mul( 0.9 );
            soldiers[i].moveSoldier( dt );
        }
    }

    void render( ){
        //printf( " rendering \n" );
        //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
        //glColor3f( 0.1f, 0.1f, 0.1f );
        glBegin( GL_LINE_LOOP );
        glVertex3f( (float)p00.x, (float)p00.y, 0 );
        glVertex3f( (float)p01.x, (float)p01.y, 0 );
        glVertex3f( (float)p11.x, (float)p11.y, 0 );
        glVertex3f( (float)p10.x, (float)p10.y, 0 );
        glEnd();

        for( int i = 0; i<nsoldiers; i++ ){
            Draw2D::drawCircle_d( soldiers[i].pos, 0.5, 8, true );
            /*
            printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n", soldiers[i].pos.x,soldiers[i].pos.y,
                           soldiers[i].vel.x,soldiers[i].vel.y,
                           soldiers[i].force.x,soldiers[i].force.y,
                           soldiers[i].willForce.x,soldiers[i].willForce.y  );
            */
        }
        //printf( "============== \n" );
        //exit(0);
    }

    Formation(){};
    Formation( int nrows_, int ncols_, SoldierType * type ){
        nrows     = nrows_;
        ncols     = ncols_;
        nsoldiers = ncols*nrows;
        nalive    = nsoldiers;
        ncapable  = nsoldiers;
        soldiers  = new Soldier[ nsoldiers ];
        setupSoldiers( type );
    }

    void setupSoldiers( SoldierType * type ){
        for( int i=0; i<nsoldiers; i++ ){
            soldiers[i].type  = type;
            soldiers[i].setMass( 1.0 );
            soldiers[i].vel      .set( 0.0 );
            soldiers[i].willForce.set( 0.0 );
            soldiers[i].force    .set( 0.0 );
        }
    }

    void deploySoldiers( ){
        double dlf = length / nrows;
        double dfw = width  / ncols;
        int i=0;
        double cfw = 0.5*dfw;
        for( int icol=0; icol<ncols; icol++ ){
            double clf = 0.5*dlf;
            for( int irow=0; irow<nrows; irow++ ){
                soldiers[i].pos.set_lincomb( 1, cfw+randf(-0.1,-0.1), clf+randf(-0.1,-0.1), p11, dirFw, dirLf );
                clf += dlf;
                i++;
            }
            cfw += dfw;
        }
/*
        for( int i=0; i<nsoldiers; i++ ){
            printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",
                soldiers[i].pos.x,soldiers[i].pos.y,
                soldiers[i].vel.x,soldiers[i].vel.y,
                soldiers[i].force.x,soldiers[i].force.y,
                soldiers[i].willForce.x,soldiers[i].willForce.y  );
        }
*/
    }
    ~Formation(){
        if( soldiers != NULL ) delete soldiers;
    }

};

#endif
