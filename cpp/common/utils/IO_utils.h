
#ifndef  IO_utils_h
#define  IO_utils_h

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <cstring>

#include <string>
#include <vector>

#include "Vec2.h"
#include "Vec3.h"


const int N_CHAR_TMP = 256;

inline int fileExist(const char * fname ){
    FILE *file;
    if ( (file = fopen(fname, "r")) ) {
        fclose(file);
        return 1;
    } else {
        return 0;
    }
}

#include <vector>
#include <unistd.h>
#include <dirent.h>

//#include "Tree.h"

// list files in directory
//  https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c

inline int listDirContaining( char * dirName, char * fname_contains, std::vector<std::string>& fnames_found ){
    DIR *dir=NULL;
    int n=0;
    struct dirent *ent=NULL;
    int i=0;
    if ( (dir = opendir( dirName )) != NULL) {
        while ( (ent = readdir (dir)) != NULL) {
            char* found = strstr( ent->d_name, fname_contains );
            if( found ){
                printf("%i %s\n", i, ent->d_name);
                fnames_found.push_back( ent->d_name );
                i++;
            }
        }
        n++;
        closedir(dir);
    } else {
        printf("Cannot open directory %s \n", dirName );
        return -1;
    }
    return i;
}

/*
int dir2tree(TreeViewTree& node, char * name, int level ){

    if (niters >100) return -1;
    niters++;

    if((name[0]=='.'))return 0;

    for(int i=0; i<level; i++) printf("_");

    node.content.caption = name;
    DIR *dir=NULL;
    struct dirent *ent=NULL;

    if( chdir(name)==0 ){
    //if( (dir = opendir( name )) != NULL){
        dir = opendir( "." );
        printf("dir '%s' | %i \n", name, level );
        while( (ent = readdir(dir)) != NULL){
            node.branches.push_back( TreeViewTree() );
            dir2tree( node.branches.back(), ent->d_name, level+1 );
        }
        closedir(dir);
        chdir("..");
    }else{
        printf("leaf '%s' | %i \n", name, level );
    }
    return 0;
}
*/

template <typename Func>
int processFileLines( char * fname, Func func ){
FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    int n=0;
    while ( fgets( str , nbuff, pFile) != NULL ){
        if (str[0]=='#') continue;
        func( str );
        n++;
    }
    fclose(pFile);
    return n;
}

inline char* stripWhite( char* s ){
    for(;;s++){
        char c=*s;
        if(c=='\0') return s;
        if(c>=33) break;
    }
    for( char* s_=s;;s_++ ){
        if(*s_<33){
            *s_='\0';
            return s;
        }
    }
}

inline int strcmp_noWhite( const char * s1, const char * s2 ){
    while( true ){
        char c=*s1;
        if(c=='\0') return -1;
        if(c>=33) break;
        s1++;
    }
    while( true ){
        char c2=*s2; if(c2=='\0') break;
        char c1=*s1;
        int d = c1-c2;
        if(d) return d;
        s1++;s2++;
    }
    return 0;
}

inline int str2enum( const char * str, int nnames, const char **names ){
    for(int i=0; i<nnames; i++ ){
        //printf( ">>%s<< >>%s<<\n", str, names[i] );
        if( strcmp_noWhite( str, names[i] ) == 0 ) return i;
    }
    return -1;
}

inline int saveBin( char *fname, int n, char * data ){
    FILE *ptr_myfile=0;
    ptr_myfile=fopen( fname,"wb");
    if (!ptr_myfile){ printf("Unable to open file! \n"); return -1; }
    int nchar = 1024;
    for( int i=1; i<=n; i+=nchar ){
        int len = nchar;
        if( (n-i)<nchar ) len = (n-i);
        fwrite( data+i, sizeof(char), len, ptr_myfile);
    }
    fclose(ptr_myfile);
    return 0;
}

inline int loadBin( char *fname, int n, char * data ){
    FILE *ptr_myfile;
    ptr_myfile=fopen( fname,"rb");
    if (!ptr_myfile){ printf("Unable to open file! \n"); return -1; }
    int nchar = 1024;
    for( int i=1; i<=n; i+=nchar ){
        int len = nchar;
        if( (n-i)<nchar ) len = (n-i);
        fread( data+i, sizeof(char), len, ptr_myfile);
    }
    fclose(ptr_myfile);
    return 0;
}

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

inline int fileGetNextKey( FILE  *fptr, char * keystr, char * tmp ){
	int nch = strlen(keystr);
    while(fgets(tmp, N_CHAR_TMP, fptr)){  if( strncmp(tmp,keystr,nch)==0 ) return ftell(fptr); };
    return -1;
}

inline char * fileCut( FILE * fptr, int ibeg, int iend ){
    int nch     = iend-ibeg;
    char * buff = new char[nch+1];
    fseek(fptr, ibeg, SEEK_SET);
    fread(buff, nch, 1, fptr);
    buff[nch]='\0';
    return buff;
}

// cut piece out of file
inline char* fileGetSection(char const  *fname, char * keystr, char * endstr ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
    char *buff = NULL;
	char      tmp[N_CHAR_TMP];
    int ibeg = fileGetNextKey(fptr, keystr, tmp);
    if( ibeg>=0 ){
        int iend = fileGetNextKey(fptr, endstr, tmp);
        if(iend>=0){           // cut (ibeg,iend)
            buff = fileCut(fptr, ibeg, iend-strlen(endstr)-1 );
        }
    }
    fclose (fptr);
	return buff;
}

inline int whichKey( char * tmp, int nkey, char ** keys ){
    for(int ikey=0; ikey<nkey; ikey++){
        char * key = keys[ikey];
        int i=0;
        bool match=true;
        while(key[i]!='\0'){ if(key[i]!=tmp[i])match=false; i++; }
        if(match) return ikey;
    }
    return -1;
}

inline char ** fileGetSections(char const  *fname, int nkey, char ** keys, char *begstr ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
	char      tmp[N_CHAR_TMP];
	int nb = strlen(begstr);
	char** result = new char*[nkey];
	for(int ikey=0; ikey<nkey; ikey++){ result[ikey] = NULL; }
	int ikey=-1,i0=-1,i1;
	while( (i1=fileGetNextKey( fptr, begstr, tmp ))>=0 ){
        //if((ikey>=0)&&(i0>=0)){
        //    //printf(" ikey %i i0  %i i1 %i \n", ikey, i0, i1 );
        //    result[ikey] = fileCut( fptr, i0, i1 );
        //}
        if((ikey>=0)&&(i0>=0)){ result[ikey] = fileCut( fptr, i0, i1 ); }
        ikey = whichKey( tmp+nb, nkey, keys );
        i0=i1;
	};
	fseek(fptr, 0, SEEK_END);
	if((ikey>=0)&&(i0>=0)) result[ikey] = fileCut( fptr, i0, ftell(fptr) ); // seaction at end of file
	return result;
}

inline int checkNullSections( int n, char ** sections ){
    for(int i=0; i<n; i++){ if(sections[i]==NULL) return i; }
    return 0;
}

inline void saveStr( char * fname, char * str ){
    int n = strlen(str);
    FILE  *fptr = fopen(fname, "wb");
    fwrite( str, n, 1, fptr);
    fclose(fptr);
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
