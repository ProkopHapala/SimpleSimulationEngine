
void setPixelsFunction( SDL_Surface* surface, int ix0, int iy0, int ix1, int iy1, Function2d func ){
	SDL_LockSurface( surface  );
	int nx = ix1-ix0;
	int ny = iy1-iy0;
	for (int iy=iy0; iy<iy1; iy++){
		Uint32 * pixel  = ((Uint32*)surface->pixels) + iy*surface->w + ix0;
		double y = i2y( iy );
		for (int ix=ix0; ix<ix1; ix++){
			double x = i2x( ix );
			double f = func( x, y );
			Uint8  fc    = ( (Uint8)(255.0*f) ) & 0xFF; 
			Uint32 color = (fc<<16)|(fc<<8)|fc | 0xFF000000;
			*pixel = color;
			pixel++;
		}
	}
	SDL_UnlockSurface( surface  );
}
