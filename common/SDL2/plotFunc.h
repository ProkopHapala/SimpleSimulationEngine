

void plotHline( SDL_Renderer* render, double y ){
	SDL_RenderDrawLine( render, 0, y2i(-y), SCREEN_WIDTH, y2i(-y) );
}

void plotVline( SDL_Renderer* render, double x ){
	SDL_RenderDrawLine( render, x2i(x), 0, x2i(x), SCREEN_HEIGHT );
}

void plotAxes( SDL_Renderer* render ){
	plotHline( render, 0 );
	plotVline( render, 0 );
}


void plotFunc( SDL_Renderer* render, int n, double * xs, double * ys, double yscale ){
	int oix=x2i(  xs[0]        );
	int oiy=y2i( -ys[0]*yscale );
	for (int i=1; i<n; i++){
		int ix = x2i(  xs[i]        );
		int iy = y2i( -ys[i]*yscale );
		SDL_RenderDrawLine( render, oix, oiy, ix, iy );
		//printf( " %i  :  %f %f %i %i \n", i, xs[i], ys[i] , ix, iy  );
		oix = ix;
		oiy = iy;
	}
}
