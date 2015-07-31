

	void render_1( Vec2f p1, Vec2f p2, Vec2f p3 ){
		glVertex3f( p1.x, p1.y, 0 ); 
		glVertex3f( p2.x, p2.y, 0 ); 
		glVertex3f( p3.x, p3.y, 0 ); 
	}

	void render_2( Vec2f p00, Vec2f p01, Vec2f p10, Vec2f p11 ){
		render_1( p00, p01, p10 );
		render_1( p11, p01, p10 );
	}

	void render_3( Vec2f p0, Vec2f p1, Vec2f p2, Vec2f p3, Vec2f p4 ){
		render_1( p0, p1, p2 );
		render_1( p0, p2, p3 );
		render_1( p0, p3, p4 );
	}

	void render_4( Vec2f p0, Vec2f p1, Vec2f pL, Vec2f pR ){
		Vec2f p0L,p1L,p0R,p1R;
		p0L.set_add( p0, pL ); p0L.mul( 0.5 );
		p0R.set_add( p0, pR ); p0R.mul( 0.5 );
		render_1( p0, p0L, p0R );
		p1L.set_add( p1, pL ); p1L.mul( 0.5 );
		p1R.set_add( p1, pR ); p1R.mul( 0.5 );
		render_1( p1, p1L, p1R );
		render_1( p0L, p0R, p1L );
		render_1( p1L, p1R, p0R );
	}

	void render_block( Vec2f p00, Vec2f p01, Vec2f p10, Vec2f p11,  int code ){
		//printf( " %i \n", code );
		switch( code ){
			case 1 : render_1 ( p00, (p00+p01)*0.5, (p00+p10)*0.5 );           break;
			case 2 : render_1 ( p01, (p01+p00)*0.5, (p01+p11)*0.5 );           break;
			case 3 : render_2 ( p00, p01, (p00+p10)*0.5, (p01+p11)*0.5 );      break;
			case 4 : render_1 ( p10, (p10+p00)*0.5, (p10+p11)*0.5 );           break;
			case 5 : render_2 ( p00, p10, (p00+p01)*0.5, (p10+p11)*0.5 );      break;
			case 6 : render_4 ( p01, p10, p00, p11                     );      break;
			case 7 : render_3 ( p00, p01, (p01+p11)*0.5, (p10+p11)*0.5, p10 ); break;
			case 8 : render_1 ( p11, (p11+p01)*0.5, (p11+p10)*0.5 );           break;
			case 9 : render_4 ( p00, p11, p01, p10 );                          break;
			case 10: render_2 ( p01, p11, (p00+p01)*0.5, (p10+p11)*0.5 );      break;
			case 11: render_3 ( p01, p00, (p00+p10)*0.5, (p11+p10)*0.5, p11 ); break;
			case 12: render_2 ( p10, p11, (p00+p10)*0.5, (p01+p11)*0.5 );      break;
			case 13: render_3 ( p10, p00, (p00+p01)*0.5, (p11+p01)*0.5, p11 ); break;
			case 14: render_3 ( p11, p01, (p01+p00)*0.5, (p10+p00)*0.5, p10 ); break; 
			case 15: render_2 ( p00, p10, p01, p11 );                          break;
		}
	}


class MapTile2D{
	public:
	int pow;
	int n;
	int n2;

	MapCell2D * cells;
	int render_list;

	MapTile2D( int pow_ ){
		pow = pow_;
		n  = 1 << pow;
		n2 = n*n;
		cells  = new MapCell2D[ n2 ];
		render_list = 0;
	}

	inline int getIndex       ( int ix, int iy ){  return (iy<<pow) + iy; };
	inline     MapCell2D * get( int ix, int iy ){  return &cells[  getIndex(ix,iy) ];  };

	int make_render( float x0, float y0, float dx, float dy, MapTile2D * tL, MapTile2D * tD, MapTile2D * tLD ){
		if( render_list > 0 ) {	glDeleteLists( render_list, 1 );	}
		render_list = glGenLists(1);
		glNewList( render_list, GL_COMPILE );
		glBegin( GL_TRIANGLES );
		glShadeModel ( GL_FLAT );
		glColor3f( 0.8, 0.6, 0.3 );
		//printf( "  %i %i \n ", pow, n );
		for( int iy=0; iy<n; iy++ ){
			int icelly = iy<<pow;
			for( int ix=0; ix<n; ix++ ){
				int icell = icelly + ix;
				//printf( "  %i %i %i \n ", ix, iy, icell );
				Vec2f p00,p01,p10,p11;
				p00.set( x0 + ix*dx, y0 + iy*dx );
				p01.set( p00.x+dx, p00.y    );
				p10.set( p00.x   , p00.y+dx );
				p11.set( p00.x+dx, p00.y+dy );
				float v00,v01,v10,v11;
				v00 = cells[ icell   ].val;
				bool xedge = (ix==n);
				bool yedge = (iy==n);
				if( xedge         ){ if( tL ==NULL ){ v01=-1; }else{ v01 = tL ->cells[ icelly ].val; } }else{ v01 = cells[ icell+1   ].val; }
				if( yedge         ){ if( tD ==NULL ){ v10=-1; }else{ v10 = tD ->cells[     ix ].val; } }else{ v10 = cells[ icell+n   ].val; }
				if( xedge&&yedge  ){ if( tLD==NULL ){ v11=-1; }else{ v11 = tLD->cells[      0 ].val; } }else{ v11 = cells[ icell+n+1 ].val; }
				int code = 0;
				if( v00 > 0 ) code = code | 0b0001;
				if( v01 > 0 ) code = code | 0b0010;
				if( v10 > 0 ) code = code | 0b0100;
				if( v11 > 0 ) code = code | 0b1000;
				render_block( p00, p01, p10, p11, code );
			}	
		}		
		//printf( " render_list %i \n ", render_list );
		glEnd();
		glBegin( GL_QUADS );
		float d = 0.05;
		//printf( "  %i %i \n ", pow, n );
		for( int iy=0; iy<n; iy++ ){
			int icelly = iy<<pow;
			for( int ix=0; ix<n; ix++ ){
				int icell = icelly + ix;
				//printf( "  %i %i %i \n ", ix, iy, icell );
				float c = (cells[ icell   ].val + 0.5)*2;
				glColor3f( c, c, c );
				Vec2f p00,p01,p10,p11;
				p00.set( x0 + ix*dx, y0 + iy*dx );
				glVertex3f( p00.x-d, p00.y-d, 1 );
				glVertex3f( p00.x+d, p00.y-d, 1 );
				glVertex3f( p00.x+d, p00.y+d, 1 );
				glVertex3f( p00.x-d, p00.y+d, 1 ); 
			}	
		}		
		glEnd();
		glEndList();
		//printf( " render_list %i \n ", render_list );
		return render_list;
	}

	int make_render_smooth( float x0, float y0, float dx, float dy, MapTile2D * tL, MapTile2D * tD, MapTile2D * tLD ){
		if( render_list > 0 ) {	glDeleteLists( render_list, 1 );	}
		render_list = glGenLists(1);
		glNewList( render_list, GL_COMPILE );
		glShadeModel ( GL_SMOOTH           );
		//printf( "  %i %i \n ", pow, n );
		for( int iy=0; iy<n; iy++ ){
			int icelly = iy<<pow;
			for( int ix=0; ix<n; ix++ ){
				int icell = icelly + ix;
				//printf( "  %i %i %i \n ", ix, iy, icell );
				Vec2f p00,p01,p10,p11,pc;
				p00.set( x0 + ix*dx, y0 + iy*dx );
				p01.set( p00.x+dx, p00.y    );
				p10.set( p00.x   , p00.y+dx );
				p11.set( p00.x+dx, p00.y+dy );
				float v00,v01,v10,v11;
				v00 = cells[ icell   ].val;
				bool xedge = (ix==n);
				bool yedge = (iy==n);
				if( xedge         ){ if( tL ==NULL ){ v01=-1; }else{ v01 = tL ->cells[ icelly ].val; } }else{ v01 = cells[ icell+1   ].val; }
				if( yedge         ){ if( tD ==NULL ){ v10=-1; }else{ v10 = tD ->cells[     ix ].val; } }else{ v10 = cells[ icell+n   ].val; }
				if( xedge&&yedge  ){ if( tLD==NULL ){ v11=-1; }else{ v11 = tLD->cells[      0 ].val; } }else{ v11 = cells[ icell+n+1 ].val; }
				pc.set( (p00+p01+p10+p11)*0.25 );
				float vc = 0.25*(v00 + v01 + v10 + v11);
				glBegin( GL_TRIANGLE_FAN );
					glColor3f( vc, vc, vc );    glVertex3f( pc.x,  pc.y,  0 ); 
					glColor3f( v00, v00, v00 ); glVertex3f( p00.x, p00.y, 0 ); 
					glColor3f( v01, v01, v01 ); glVertex3f( p01.x, p01.y, 0 ); 
					glColor3f( v11, v11, v11 ); glVertex3f( p11.x, p11.y, 0 ); 
					glColor3f( v10, v10, v10 ); glVertex3f( p10.x, p10.y, 0 ); 
					glColor3f( v00, v00, v00 ); glVertex3f( p00.x, p00.y, 0 ); 
				glEnd();
			}	
		}		
		//printf( " render_list %i \n ", render_list );
		glBegin( GL_QUADS );
		float d = 0.05;
		//printf( "  %i %i \n ", pow, n );
		for( int iy=0; iy<n; iy++ ){
			int icelly = iy<<pow;
			for( int ix=0; ix<n; ix++ ){
				int icell = icelly + ix;
				//printf( "  %i %i %i \n ", ix, iy, icell );
				float c = cells[ icell   ].val;
				glColor3f( c, c, c );
				Vec2f p00,p01,p10,p11;
				p00.set( x0 + ix*dx, y0 + iy*dx );
				glVertex3f( p00.x-d, p00.y-d, 1 );
				glVertex3f( p00.x+d, p00.y-d, 1 );
				glVertex3f( p00.x+d, p00.y+d, 1 );
				glVertex3f( p00.x-d, p00.y+d, 1 ); 
			}	
		}		
		glEnd();
		glEndList();
		//printf( " render_list %i \n ", render_list );
		return render_list;
	}

};
