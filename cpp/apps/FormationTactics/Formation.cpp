
#include "Formation.h" // THE HEADER

#include <algorithm>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

void Formation::update_bbox( ){
    Vec2d p = soldiers[0].pos;
    cog.set( p );
    bbox.a.set(p); bbox.b.set(p);
    for(int i=1; i<nCapable; i++){
        p = soldiers[i].pos;
        cog.add(p);
        if( p.x < bbox.x0 ){ bbox.x0 = p.x; }
        if( p.y < bbox.y0 ){ bbox.y0 = p.y; }
        if( p.x > bbox.x1 ){ bbox.x1 = p.x; }
        if( p.y > bbox.y1 ){ bbox.y1 = p.y; }
    }
    cog.mul( 1.0/nCapable );
    bbox.a.add(-bboxMargin); bbox.b.add(bboxMargin);
}

void Formation::moveBy( const Vec2d& dpos ){
    center.add( dpos );
    setCenterRot( center, dirFw );
}

int Formation::interact( Formation * fb ){
    //if ( fb == NULL ) return;
    if ( bbox.notOverlaps( fb->bbox ) ) return 0;
    bool enemy = ( fb->faction != faction );
    for( int i=0; i<fb->nCapable; i++ ){
        Soldier * si = fb->soldiers + i;
        for( int j=0; j<nCapable; j++ ){
            Soldier * sj = soldiers + j;
            //si->enemy_interaction( sj, meele&&enemy );
            if( enemy ) {
                si->enemy_interaction ( sj, melee );
                //si->enemy_interaction ( sj, false );
            }else{
                si->friend_interaction( sj );
            }

        }
    }
    return nCapable * fb->nCapable;
}

int Formation::interactInside( ){
    for( int i=0; i<nCapable; i++ ){
        Soldier * si = soldiers + i;
        //for( int j=0; j<i; j++ ){
        for( int j=0; j<nCapable; j++ ){
            Soldier * sj = soldiers + j;
            si->friend_interaction( sj );
        }
    }
    return nCapable*nCapable;
}

void Formation::setTarget( const Vec2d& target ){
    //printf( " setTarget (%3.3f,%3.3f) \n", target.x, target.y );
    Vec2d d;
    d.set_sub( target, center ); d.normalize();
    //dirFw.set( d ); dirFw.normalize(); dirLf.set_perp( dirFw );
    setCenterRot( center, d );
    p00target.set_add_mul( target, dirLf,  length * 0.5d );
    p01target.set_add_mul( target, dirLf, -length * 0.5d );
    movingToTarget = true;
}

void Formation::jumpToTarget( ){
    setEnds( p00target, p01target, width );
}

void Formation::leaveMenBehind( ){
    for( int i=0; i<nCapable; i++ ){
        Vec2d d; d.set_sub( soldiers[i].pos, center );
        double llf = -dirLf.dot( d );
        double lfw = -dirFw.dot( d );
        if( ( fabs(llf)>(length+bboxMargin) ) || ( fabs(lfw)>(width+bboxMargin) ) ){
            //printf( "soldier %i abandoned \n" );
            soldiers[i].impair_mask |= 8;
        }
    }
}

bool Formation::eliminateInvalids( ){
    int j = nCapable-1;
    bool change = false;
    for( int i=0; i<nCapable; i++ ){
        if( soldiers[i].impair_mask >=8 ){
            change = true;
            while( soldiers[j].impair_mask >=8 ){ j--; }
            if( j>i ){
                std::swap( soldiers[i], soldiers[j] );
            }else{
                j = i;
                break;
            }
        }
    }
    if( change ){
        //printf( " nCapable %i -> %i \n", nCapable, j );
        //for( int i=0; i<nSoldiers; i++ ){ printf( " soldier %i : %i \n", i, soldiers[i].impair_mask ); } // just debug
        nCapable = j;
    }
    return change;
}

void Formation::moveToTarget( ){

    if( checkMenBehind( ) ){
        //printf( " formation %i cannot move, men stuck ! \n", id );
        if( shouldLeaveMenBehind ){
            //printf( " => leaving men behind ! \n", id );
            leaveMenBehind( );
        }
        return;
    }

    checkTarget( );
    if( movingToTarget ){
        double speed = 0.01;
        Vec2d d1,d0;
        d0.set_sub( p00target, p00 ); d0.normalize(); p00.add_mul( d0, speed );
        d1.set_sub( p01target, p01 ); d1.normalize(); p01.add_mul( d1, speed );
        setEnds( p00, p01, width );
    }
}

void Formation::applyWillForce( Soldier& soldier ){
    //if( soldier.impair_mask < 8 ){
        Vec2d d;
        d.set_sub( soldier.pos, p00 );
        double llf = -dirLf.dot( d );
        if      ( llf < 0      ){ soldier.willForce.add_mul( dirLf, -willSaturation( -llf        *kLength)); }
        else if ( llf > length ){ soldier.willForce.add_mul( dirLf,  willSaturation( (llf-length)*kLength)); }
        double lfw = -dirFw.dot( d );
        if      ( lfw < 0      ){ soldier.willForce.add_mul( dirFw, -willSaturation( -lfw        *kWidth )); }
        else if ( lfw > width  ){ soldier.willForce.add_mul( dirFw,  willSaturation( (lfw-width )*kWidth )); }
        //printf( "- %3.3f %3.3f   %3.3f %3.3f \n", llf, lfw, length, width  );
    //}
}

void Formation::applyWillForce( ){
    for( int i = 0; i<nCapable; i++ ){
        applyWillForce( soldiers[i] );
        //printf( " (%3.3f,%3.3f)\n", i, soldiers[i].willForce.x, soldiers[i].willForce.y );
    }
}

void Formation::clean_temp(){
    for( int i = 0; i<nCapable; i++ ){
        soldiers[i].clean_temp();
    }
}

void Formation::setEnds( const Vec2d& pL_, const Vec2d& pR_, double width_ ){
    width = width_;
    p00.set( pL_ ); p01.set( pR_ );
    center.set_lincomb( 0.5d, p00, 0.5d, p01 );
    Vec2d d; d.set_sub( p00, p01 );
    length = d.norm();
    dirLf.set_mul ( d, 1/length  );
    dirFw.set_perp( dirLf );  dirFw.mul(-1);
    p10.set_add_mul( p00, dirFw, -width );
    p11.set_add_mul( p01, dirFw, -width );
    //printf( " p00 (%3.3f,%3.3f) p01 (%3.3f,%3.3f) dirLf (%3.3f,%3.3f) dirFw (%3.3f,%3.3f) \n", p00.x,p00.y, p01.x,p01.y, dirLf.x,dirLf.y, dirFw.x,dirFw.y );
    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
}


void Formation::setCenterRot( const Vec2d& center_, const Vec2d& dirFw_ ){
    center.set( center_ );
    dirFw .set( dirFw_  );
    dirLf .set_perp ( dirFw );
    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) %3.3f %3.3f \n",center.x,center.y, dirFw.x,dirFw.y, width, length );
    p00.set_add_mul( center, dirLf,  length * 0.5d );
    p01.set_add_mul( center, dirLf, -length * 0.5d );
    p10.set_add_mul( p00,    dirFw, -width );
    p11.set_add_mul( p01,    dirFw, -width );
    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
    //exit(0);
}

void Formation::update( double dt ){
    moveToTarget( );
    for( int i = 0; i<nCapable; i++ ){
        //soldiers[i].rot.add_mul( dirFw, 0.1 );
        //soldiers[i].rot.normalize();
        //soldiers[i].vel.mul( 0.9 );
        applyWillForce( soldiers[i] );
        soldiers[i].attentionDir.add_mul( dirFw, 0.1 );
        soldiers[i].update      ( dt, tAttack );
    }
}

void Formation::render( const Vec3f& color, int view_type ){
    //printf( " rendering \n" );
    //printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) \n",p00.x,p00.y, p01.x,p01.y, p10.x,p10.y, p11.x,p11.y );
    //glColor3f( faction->color.x, faction->color.y, faction->color.z );

    glColor3f( color.x, color.y, color.z );

    glBegin( GL_LINE_LOOP );
    glVertex3f( (float)p00.x, (float)p00.y, 0 );
    glVertex3f( (float)p01.x, (float)p01.y, 0 );
    glVertex3f( (float)p11.x, (float)p11.y, 0 );
    glVertex3f( (float)p10.x, (float)p10.y, 0 );
    glEnd();

    Draw2D::drawRectangle_d   ( bbox.a, bbox.b, false );
    Draw2D::drawPointCross_d( cog, 0.5 );

    if( movingToTarget ){
        Draw2D::drawLine_d( p00target, p01target );
        Draw2D::drawLine_d( center, (p01target+p00target)*0.5 );
    }

    for( int i = 0; i<nCapable; i++ ){
        Soldier& si = soldiers[i];
        float c;
        Vec3f cv;
        switch( view_type ){
            case VIEW_INJURY:  cv = si.impair2color(); glColor3f( cv.x, cv.y, cv.z ); break;
            case VIEW_STAMINA: c  = si.stamina;               glColor3f( c,c,c ); break;
            case VIEW_CHARGE:  c  = si.time_buf/5.0; if(c>0){ glColor3f( 1-c,1-c,1 ); }else{ glColor3f( 1,1+c,1+c); } break;
            case VIEW_MORAL:   c  = si.moral;                 glColor3f( c,c,c ); break;
        };
        //Draw2D::drawCircle_d( soldiers[i].pos, 0.5, 8, true );
        Draw2D::drawCircle_d( soldiers[i].pos, 0.25, 8, true );
        //Draw2D::drawLine_d  ( soldiers[i].pos, soldiers[i].pos );
        Draw2D::drawVecInPos_d( soldiers[i].rot*soldiers[i].type->melee_range, soldiers[i].pos );



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

char* Formation::reportStatus( char * sout ){
    sout += sprintf( sout, "able %i alive %i of %i \n", nCapable,nAlive,nSoldiers );
    sout += sprintf( sout, "order             %f\n", order                  );
    sout += sprintf( sout, "moral             %f\n", moral                  );
    sout += sprintf( sout, "stamina           %f\n", stamina                );
    sout += sprintf( sout, "movingToTarget    %c\n", movingToTarget?'T':'F' );
    //sout += sprintf( sout, "stamina_regain   %f\n", stamina_regain         );
    return sout;
}

char* Formation::reportSetup( char * sout ){
    sout += sprintf( sout, "Formation %i  %s \n", id, name.c_str()  );
    sout += sprintf( sout, "n(row,col)    %i(%i,%i) \n", nSoldiers, nrows, ncols );
    sout += sprintf( sout, "(lengh,width)   (%f,%f)  \n", length,width );
    sout += sprintf( sout, "(klengh,kwidth) (%f,%f)  \n", kLength, kWidth );
    sout += sprintf( sout, "maxWill        %f\n", maxWill           );
    sout += sprintf( sout, "bboxMargin     %f\n", bboxMargin        );
    sout += sprintf( sout, "maxBbox2       %f\n", maxBbox2          );
    sout += sprintf( sout, "tAttack        %f\n", tAttack           );
    sout += sprintf( sout, "melee          %c\n", melee?'T':'F'     );
    sout += sprintf( sout, "LeaveMenBehind %c\n", shouldLeaveMenBehind?'T':'F' );
    return sout;
};


Formation::Formation( int id_, int nrows_, int ncols_, SoldierType * type, Faction * faction_ ){
    id = id_;
    //name      = name_;
    faction   = faction_;
    nrows     = nrows_;
    ncols     = ncols_;
    nSoldiers = ncols*nrows;
    nAlive    = nSoldiers;
    nCapable  = nSoldiers;
    soldiers  = new Soldier[ nSoldiers ];
    setupSoldiers( type );
}

void Formation::setupSoldiers( SoldierType * type ){
    bboxMargin = type->melee_range;
    for( int i=0; i<nSoldiers; i++ ){
        soldiers[i].type  = type;
        soldiers[i].setMass( type->mass );
        soldiers[i].vel      .set( 0.0 );
        soldiers[i].willForce.set( 0.0 );
        soldiers[i].force    .set( 0.0 );
        soldiers[i].attentionDir.set( 0.0 );
    }
}

void Formation::deploySoldiers( ){
    //printf( "length %f width %f \n", length, width  );
    double dlf = length / ncols;
    double dfw = width  / nrows;
    int i=0;
    double cfw = 0.5*dfw;
    for( int irow=0; irow<nrows; irow++ ){
        double clf = 0.5*dlf + ( (irow&1) -0.5d );
        for( int icol=0; icol<ncols; icol++ ){
            //soldiers[i].pos.set_lincomb( 1, cfw+randf(-0.1,-0.1), clf+randf(-0.1,-0.1), p11, dirFw, dirLf );
            soldiers[i].pos.set_lincomb( 1, cfw, clf, p11, dirFw, dirLf );
            clf += dlf;
            i++;
        }
        cfw += dfw ;
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
