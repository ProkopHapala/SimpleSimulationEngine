
#ifndef  MusicUtils_h
#define  MusicUtils_h

#include "SDL2/SDL_mixer.h"
#include "Fourier.h"

namespace Music{

class Spectrum{ public:

    int nwave=0;
    int nhist=0;

    double*    wave=0; // complex numbers

    double*    hist   =0;
    double*    histOld=0;
    double*    histVel=0;


    int position=0;
    bool need_refresh=0;


    int    done=0, bits=0, which=0, sample_size=0, rate=0;
    int    audio_rate,audio_channels;
    Uint16 audio_format;
    //float dy;

    void realloc( int nwave_, int nhist_ ){
        nwave=nwave_;
        nhist=nhist_;
        _realloc( wave,    nwave );
        _realloc( hist,    nwave );
        _realloc( histOld, nwave );
        _realloc( histVel, nwave );
    }

    void clearHist(){
        for(int i=0; i<nhist; i++){ hist[i]=0; histOld[i]=0; histVel[i]=0; }
    }

    void powerSpectrum( ){
        for(int i=0; i<nhist; i++){ hist[i] = 0; }
        double dx = 2*nhist/nwave;
        for(int i=0; i<nwave/2; i++){
            int i1=i*2;
            int i2=2*nwave-i1;
            double yr1=wave[i1  ];
            double yi1=wave[i1+1];
            double yr2=wave[i2  ];
            double yi2=wave[i2+1];
            double p1 = yr1*yr1 + yi1*yi1;
            double p2 = yr2*yr2 + yi2*yi2;
            int j = (int)(i*dx);
            hist[j] += p1+p2;
        }
        for(int i=0; i<nhist; i++){ hist[i] = sqrt( hist[i] ); }
    }

    void spectrumHistDynamics( double da ){
        for(int i=0; i<nhist; i++){
            double y = hist[i];
            //printf( "%g %g \n", y, old[i] );
            if(y<histOld[i]){ histVel[i]-=da; histOld[i]+=histVel[i]; }
            if(y>histOld[i]){ histOld[i]=y;   histVel[i]=0;           }
        }
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
        printf("Opened audio at %d Hz %d bit %s, %d bytes audio buffer\n", audio_rate, bits, audio_channels>1?"stereo":"mono", nwave);

        //dy=HEIGHT/2.0/(float)(0x1<<bits);
    }


};


static void postmix_Spectrum( void *spec_, Uint8 *_stream, int _len){

    Spectrum& spec = *(Spectrum*) spec_;

	spec.position+=_len/spec.sample_size;
	// fprintf(stderr,"pos=%7.2f seconds \r",position/(float)rate);
	//printf( "position %i len %i sample_size %i \n", spec.position, _len, spec.sample_size );
	//return;
	if(spec.need_refresh) return;
	// save the stream buffer and indicate that we need a redraw
	//len=_len;

    /// ??? QUESTION ???: which-fliping is there because of stereo ?????
	//memcpy(stream[(which+1)%2],_stream,len>WIDTH*4?WIDTH*4:len);
    //spec.which=(spec.which+1)%2;

    printf( "position %i len %i sample_size %i \n", spec.position, _len, spec.sample_size );

	Sint16* stream = (Sint16*)_stream;
	//for(int i=0; i<spec.nwave; i++){

	for(int i=0; i<_len; i++){
        spec.wave[i] = (double)stream[i];
	}

	spec.need_refresh=true;

}





};

#endif

