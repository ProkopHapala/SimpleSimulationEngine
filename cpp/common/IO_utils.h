
#ifndef  IO_utils_h
#define  IO_utils_h

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <cstring>

#include "Vec2.h"
#include "Vec3.h"


inline char * fgets_comment( char * line, int num, FILE * stream ){
    constexpr int NMaxComment = 10;
    for(int i=0; i<NMaxComment; i++){
        char *str = fgets( line, num, stream );
        printf(">> %s", line );
        if( str[0] != '#' ) return str;
    }
    return NULL;
}

// A simple function that will read a file into an allocated char pointer buffer
inline  char* filetobuf(char const  *fname){
	FILE *fptr;
	long length;
	char *buf;

	fptr = fopen(fname, "rb");			// Open file for reading
	if(fptr==NULL){
        printf("Failed to load %s \n", fname );
	    return NULL;
    }
	fseek(fptr, 0, SEEK_END); 			// Seek to the end of the file
	length = ftell(fptr); 				// Find out how many bytes into the file we are
	//buf = (char*)malloc(length+1); 		// Allocate a buffer for the entire length of the file and a null terminator
	buf = new char[length+1];
	fseek(fptr, 0, SEEK_SET); 			// Go back to the beginning of the file
	fread(buf, length, 1, fptr); 		// Read the contents of the file in to the buffer
	fclose(fptr); 						// Close the file
	buf[length] = 0; 					// Null terminator
	return buf; 						// Return the buffer
}

/*
int loadColums(char const  *fname, char const  *format, ... ){
    FILE * pFile = fopen (fname,"r");
    char buff[1024];
    char * line;
    while( line = fgets( buff, 1024, pFile ) ){
        sscanf( buff, format, ... );
        printf()
    }
    va_list args;
    va_start(args, fmt);
    scanf( format );
    va_end(args);
    fclose(pFile);
}
*/

inline  int allocateIOBuffs( int nitems, char const *format, void **buffs ){
    int nbuffs = 0;
    //int ibuff  = 0;
    while (*format != '\0') {
        printf( "format %c nbuffs %i \n", *format, nbuffs );
        switch( *format ){
            case 'i': buffs[nbuffs] = new int   [nitems]; nbuffs++; break;
            case 'f': buffs[nbuffs] = new float [nitems]; nbuffs++; break;
            case 'd': buffs[nbuffs] = new double[nitems]; nbuffs++; break;
            case '2': buffs[nbuffs] = new Vec3d [nitems]; nbuffs++; break;  // TODO
            case '3': buffs[nbuffs] = new Vec2d [nitems]; nbuffs++; break;
        }
        format++;
    }
    return nbuffs;
}

inline  int loadColumns( char const  *fname, char const *format, void **buffs ){
    FILE * pFile = fopen(fname,"r");
    if( pFile == NULL ){
        printf("cannot find %s\n", fname );
        return -1;
    }
    char buff[1024];
    char * line;
    int nl;
    line = fgets( buff, 1024, pFile );
    if(line==NULL){
        printf("read nl line NULL \n");
        return -1;
    }
    sscanf( buff,"%i", &nl );
    printf(" nl = %i \n", nl);
    allocateIOBuffs( nl, format, buffs );
    for(int il=0; il<nl; il++){
        line = fgets( buff, 1024, pFile );
        //while( line = fgets( buff, 1024, pFile ) ){
        //printf("%s \n", line );
        int ib = 0;
        const char *formati = format;
        char *tok = strtok(line, " \t");
        while (tok != NULL) {
            //my_array[i++] = atof(tok);
            //printf( "%s   %c \n", tok, *formati );
            switch( *formati ){
                case 'i': ((int   *)buffs[ib])[il]=atoi(tok); ib++; break;
                case 'f': ((float *)buffs[ib])[il]=atof(tok); ib++; break;
                case 'd': ((double*)buffs[ib])[il]=atof(tok); ib++; break;
                case '2':{
                    Vec2d& v2 = ((Vec2d*)buffs[ib])[il];
                    v2.x = atof(tok); tok = strtok(NULL, " \t");
                    v2.y = atof(tok);
                    }; ib++; break;
                case '3':{
                    //printf(" Vdfsdfsdf455464 \n");
                    Vec3d& v3 = ((Vec3d*)buffs[ib])[il];
                    v3.x = atof(tok); tok = strtok(NULL, " \t"); // printf("tok_y %s \n",tok);
                    v3.y = atof(tok); tok = strtok(NULL, " \t"); // printf("tok_z %s \n",tok);
                    v3.z = atof(tok);
                    }; ib++; break;
            }
            tok = strtok(NULL, " \t");
            formati++;
        }
    }
    return nl;
}


#endif
