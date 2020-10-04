#include "SDL2/SDL.h"
#include "SDL2/SDL_mixer.h"

#include "Draw2D.h"
#include "AppSDL2OGL.h"


#include "Fourier.h"



/*

use this function to vizualize music data

void Mix_SetPostMix();

https://www.libsdl.org/projects/SDL_mixer/docs/SDL_mixer_frame.html



*/



//static const char *MY_COOL_MP3 = "cool_tunes.mp3";
static const char *music_file_name = "common_resources/02-Lazenca-SaveUs.mp3";
//static const char *music_file_name = "Yanni - Reflections Of Passion.mp3";



// set this to any of 512,1024,2048,4096
// the lower it is, the more FPS shown and CPU needed
#define BUFFER 1024
#define WIDTH  640 // NEVER make this be less than BUFFER!
#define HEIGHT 480
#define HEIGHT2 (HEIGHT/2)
#define HEIGHT4 (HEIGHT/4)
#define Y(sample) (((sample)*HEIGHT)/4/0x7fff)

char   post_state_buff[1024];
Sint16 stream[2][BUFFER*2*2];

double fft_buff[BUFFER*2+32]; // complex numbers
const int nbinHist  = 32;
double    powerHist   [nbinHist];
double    powerHistTop[nbinHist];
double    powerHistVel[nbinHist];


int len=BUFFER*2*2, done=0, need_refresh=0, bits=0, which=0, sample_size=0, position=0, rate=0;

int    audio_rate,audio_channels;
Uint16 audio_format;
float dy;

static void postmix( void *udata, Uint8 *_stream, int _len){

	position+=_len/sample_size;
	// fprintf(stderr,"pos=%7.2f seconds \r",position/(float)rate);
	printf( "position %i len %i sample_size %i \n", position, _len, sample_size );
	//return;
	if(need_refresh) return;
	// save the stream buffer and indicate that we need a redraw
	len=_len;
	// ??? QUESTION ???: which-fliping is there because of stereo ?????
	memcpy(stream[(which+1)%2],_stream,len>WIDTH*4?WIDTH*4:len);
	which=(which+1)%2;
	need_refresh=1;
}


void mixer_info(){
	// print out some info on the formats this run of SDL_mixer supports */
	{
		int i,n=Mix_GetNumChunkDecoders();
		printf("There are %d available chunk(sample) decoders:\n",n);
		for(i=0; i<n; ++i) printf("	%s\n", Mix_GetChunkDecoder(i));
		n = Mix_GetNumMusicDecoders();
		printf("There are %d available music decoders:\n",n);
		for(i=0; i<n; ++i) printf("	%s\n", Mix_GetMusicDecoder(i));
	}
	// print out some info on the audio device and stream */
	Mix_QuerySpec(&audio_rate, &audio_format, &audio_channels);
	bits        = audio_format&0xFF;
	sample_size = bits/8+audio_channels;
	rate        = audio_rate;
	printf("Opened audio at %d Hz %d bit %s, %d bytes audio buffer\n", audio_rate, bits, audio_channels>1?"stereo":"mono", BUFFER );

	dy=HEIGHT/2.0/(float)(0x1<<bits);
}



void spectrumHist( double nbuf, double* buf, int nbin, double* out){
    for(int i=0; i<nbin; i++){ out[i] = 0; }
    double dx = 2*nbin/nbuf;
    for(int i=0; i<nbuf/2; i++){
        int i1=i*2;
        int i2=2*nbuf-i1;
        double yr1=buf[i1  ];
        double yi1=buf[i1+1];
        double yr2=buf[i2  ];
        double yi2=buf[i2+1];
        double p1 = yr1*yr1 + yi1*yi1;
        double p2 = yr2*yr2 + yi2*yi2;
        int j = (int)(i*dx);
        out[j] += p1+p2;
    }
    for(int i=0; i<nbin; i++){ out[i] = sqrt( out[i] ); }
}

void spectrumHistDynamics( double da, int nbin, double* now, double* old, double* vel ){
    for(int i=0; i<nbin; i++){
        double y = now[i];
        //printf( "%g %g \n", y, old[i] );
        if(y<old[i]){ vel[i]-=da; old[i]+=vel[i]; }
        if(y>old[i]){ old[i]=y;   vel[i]=0;        }
    }
}



void redraw_waveform(){
	int     x;
	Sint16 *buf;
	buf         = stream[which];
	need_refresh= 0;
	// draw the wav from the saved stream buffer
	//for(x=0;x<WIDTH*2;x++){

	/*
	for(x=0;x<BUFFER;x++){
        const int  X=x>>1;
		int y = Y(buf[x]);
		glColor3f(1.,0.,0.); Draw2D::drawLine( {X,0}, {X,y} );
		//glColor3f(0.,0.,1.); Draw2D::drawLine( {X,0}, {X,h1} );
	}
	*/

	for(int i=0;i<BUFFER+32;i++){
        fft_buff[2*i  ]  = Y(buf[i]);  // Re[i]
        //fft_buff[2*i  ] = sin(2*200*M_PI*i/((double)BUFFER))*0 + sin(2*500*M_PI*i/((double)BUFFER)) + 1;
        //fft_buff[2*i  ] = 1;
        fft_buff[2*i+1] = 0;          // Im[i]
	}

	/*
	glColor3f(1.,1.,1.); Draw2D::drawLine( {0,0}, {BUFFER,0} );
	glBegin(GL_LINE_STRIP);
	glColor3f(0.,0.,1.);
    for(int i=0;i<BUFFER;i++){
        double y = fft_buff[2*i  ];
        //Draw2D::drawLine( {X,0}, {X,y} );
        glVertex2f(i,y);
	}
	glEnd();

    glBegin(GL_LINE_STRIP);
    glColor3f(1.,0.,0.);
    for(int i=0;i<BUFFER;i++){
        double y = fft_buff[2*i+1];
		//Draw2D::drawLine( {X,0}, {X,y} );
		glVertex2f(i,y);
	}
	glEnd();
    */

	FFT( fft_buff, BUFFER, 1 );

    /*
	float sc = 0.2;
	for(int i=0;i<BUFFER;i++){
        double yr = fft_buff[2*i  ];
        double yi = fft_buff[2*i+1];
		//glColor3f(1.,0.,0.); Draw2D::drawLine( {i    ,0}, {i    ,yr*sc} );
		//glColor3f(0.,0.,1.); Draw2D::drawLine( {i+0.5,0}, {i+0.5,yi*sc} );
		glColor3f(0.,1.,0.); Draw2D::drawLine( {i+0.5,0}, {i+0.5,sqrt(yr*yr + yi*yi)*sc} );
	}
    */

	spectrumHist( BUFFER, fft_buff, nbinHist, powerHist );
	spectrumHistDynamics( 5.1, nbinHist, powerHist, powerHistTop, powerHistVel );

	{float sc = 0.002;
	for(int i=0;i<nbinHist;i++){
        double yr = powerHist[i];
		glColor3f(1.,0.,0.); Draw2D::drawLine( {i    ,0}, {i    ,yr*sc} );
	}
	glBegin(GL_LINE_STRIP);
    glColor3f(0.,0.,1.);
    for(int i=0;i<nbinHist;i++){
        double y = powerHistTop[i];
		glVertex2f(i,y*sc);
	}
	glEnd();
	}

}

class TestAppCityGen : public AppSDL2OGL { public:
	//int viewList;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();

	TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCityGen::TestAppCityGen( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    int result = 0;
    int flags = MIX_INIT_MP3;

    Mix_OpenAudio(22050, AUDIO_S16SYS, 2, 640);

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        printf("Failed to init SDL\n");
        exit(1);
    }

    if (flags != (result = Mix_Init(flags))) {
        printf("Could not initialize mixer (result: %d).\n", result);
        printf("Mix_Init: %s\n", Mix_GetError());
        exit(1);
    }

    //Mix_OpenAudio(22050, AUDIO_S16SYS, 2, 640);
    Mix_Music *music = Mix_LoadMUS( music_file_name );

    mixer_info();

    for(int i=0; i<nbinHist; i++){ powerHist[i]=0; powerHistTop[i]=0; powerHistVel[i]=0; }

    camStep = 2.0;



    Mix_SetPostMix( postmix, post_state_buff );
    Mix_PlayMusic(music, 1);

    //while (!SDL_QuitRequested()) { SDL_Delay(250); }
    //Mix_FreeMusic(music);
}

void TestAppCityGen::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	redraw_waveform();

};

void TestAppCityGen::drawHUD(){}

TestAppCityGen * testApp;

int main(int argc, char **argv) {
  	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppCityGen( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;

}
