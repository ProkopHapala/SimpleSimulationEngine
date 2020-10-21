
#ifndef  IO_utils_h
#define  IO_utils_h

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <cstring>

#include <string>
#include <vector>
#include <unordered_map>

#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

const int N_CHAR_TMP = 256;

template <typename T> void toBuff  (  const T& v, void* buff, int& i){ (*(T*)(buff+i))=v; i+=sizeof(T); };
template <typename T> void fromBuff(        T& v, void* buff, int& i){ v=*((T*)(buff+i)); i+=sizeof(T); };

inline int print(       char*   v){ return printf( "%s", v ); };
inline int print(       float   v){ return printf( "%g", v ); };
inline int print(       double  v){ return printf( "%g", v ); };
inline int print(       int     v){ return printf( "%i", v ); };
// Moved to Vec2.h, Vec3.h, quaternion.h
//inline int print( const Vec2f&  v){ return printf( "%lg %g", v.x, v.y ); };
//inline int print( const Vec2d&  v){ return printf( "%lg %g", v.x, v.y ); };
//inline int print( const Vec2i&  v){ return printf( "%li %i", v.x, v.y ); };
//inline int print( const Vec3f&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
//inline int print( const Vec3d&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
//inline int print( const Vec3i&  v){ return printf( "%i %i %i", v.x, v.y, v.z ); };
//inline int print( const Quat4f& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
//inline int print( const Quat4d& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
//inline int print( const Quat4i& v){ return printf( "%i %i %i %i", v.x, v.y, v.z, v.w ); };

inline int fprint(FILE* sbuff,       char*   v){ return fprintf(sbuff, "%s ", v ); };
inline int fprint(FILE* sbuff,       float   v){ return fprintf(sbuff, "%g ", v ); };
inline int fprint(FILE* sbuff,       double  v){ return fprintf(sbuff, "%g ", v ); };
inline int fprint(FILE* sbuff,       int     v){ return fprintf(sbuff, "%i ", v ); };
inline int fprint(FILE* sbuff, const Vec2f&  v){ return fprintf(sbuff, "%lg %g ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec2d&  v){ return fprintf(sbuff, "%lg %g ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec2i&  v){ return fprintf(sbuff, "%li %i ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec3f&  v){ return fprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Vec3d&  v){ return fprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Vec3i&  v){ return fprintf(sbuff, "%i %i %i ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Quat4f& v){ return fprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int fprint(FILE* sbuff, const Quat4d& v){ return fprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int fprint(FILE* sbuff, const Quat4i& v){ return fprintf(sbuff, "%i %i %i %i ", v.x, v.y, v.z, v.w ); };

inline int sprint(char* sbuff,       char*   v){ return sprintf(sbuff, "%s ", v ); };
inline int sprint(char* sbuff,       float   v){ return sprintf(sbuff, "%g ", v ); };
inline int sprint(char* sbuff,       double  v){ return sprintf(sbuff, "%g ", v ); };
inline int sprint(char* sbuff,       int     v){ return sprintf(sbuff, "%i ", v ); };
inline int sprint(char* sbuff, const Vec2f&  v){ return sprintf(sbuff, "%lg %g ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec2d&  v){ return sprintf(sbuff, "%lg %g ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec2i&  v){ return sprintf(sbuff, "%li %i ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec3f&  v){ return sprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Vec3d&  v){ return sprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Vec3i&  v){ return sprintf(sbuff, "%i %i %i ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Quat4f& v){ return sprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int sprint(char* sbuff, const Quat4d& v){ return sprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int sprint(char* sbuff, const Quat4i& v){ return sprintf(sbuff, "%i %i %i %i ", v.x, v.y, v.z, v.w ); };

inline int sscan(char* sbuff, char*&  v){ int n; sscanf(sbuff, "%s%n", &v, &n );                                  return n; };
inline int sscan(char* sbuff, float&  v){ int n; sscanf(sbuff, "%f%n", &v , &n);                                  return n; };
inline int sscan(char* sbuff, double& v){ int n; sscanf(sbuff, "%lf%n", &v , &n);                                 return n; };
inline int sscan(char* sbuff, int&    v){ int n; sscanf(sbuff, "%i%n", &v , &n);                                  return n; };
inline int sscan(char* sbuff, Vec2f&  v){ int n; sscanf(sbuff, "%f %f%n", &v.x, &v.y , &n);                       return n; };
inline int sscan(char* sbuff, Vec2d&  v){ int n; sscanf(sbuff, "%lf %lf%n", &v.x, &v.y , &n);                     return n; };
inline int sscan(char* sbuff, Vec2i&  v){ int n; sscanf(sbuff, "%i %i%n", &v.x, &v.y , &n);                       return n; };
inline int sscan(char* sbuff, Vec3f&  v){ int n; sscanf(sbuff, "%f %f %f%n",    &v.x, &v.y, &v.z , &n);           return n; };
inline int sscan(char* sbuff, Vec3d&  v){ int n; sscanf(sbuff, "%lf %lf %lf%n",    &v.x, &v.y, &v.z , &n);        return n; };
inline int sscan(char* sbuff, Vec3i&  v){ int n; sscanf(sbuff, "%i %i %i%n",    &v.x, &v.y, &v.z , &n);           return n; };
inline int sscan(char* sbuff, Quat4f& v){ int n; sscanf(sbuff, "%f %f %f %f%n", &v.x, &v.y, &v.z, &v.w , &n);     return n; };
inline int sscan(char* sbuff, Quat4d& v){ int n; sscanf(sbuff, "%lf %lf %lf %lf%n", &v.x, &v.y, &v.z, &v.w , &n); return n; };
inline int sscan(char* sbuff, Quat4i& v){ int n; sscanf(sbuff, "%i %i %i %i%n",    &v.x, &v.y, &v.z, &v.w , &n);  return n; };

//#define OPT(OP,) OP(tok)
//#define OPT(OP,tok) OP(tok)
//#define DO_10(OP,tok,...) OP(tok) DO_9(__VA_ARGS__)

#define _sprintf( sbuff, format, ... ){ sbuff+=sprintf( sbuff, format, __VA_ARGS__); }
#define _sprint ( sbuff, tok ){ sbuff+=sprint( sbuff, tok ); }
#define _sscan  ( sbuff, tok ){ sbuff+=sscan( sbuff, tok ); }

#define _lprint(         tok        ) {         printf(        " %s: ",#tok);          print(      tok); }
#define _lfprint( fout,  tok        ) {         fprintf( fout,  " %s: ",#tok);        fprint(fout, tok); }
#define _lsprint( sbuff, tok        ) { sbuff+= sprintf( sbuff, " %s: ",#tok); sbuff+=sprint(sbuff,tok); }
#define _Lprint(         tok, s     ) {         printf(         " %s: ",s);          print(      tok); }
#define _Lfprint( fout,  tok, s     ) {         fprintf( fout,  " %s: ",s);        fprint(fout, tok); }
#define _Lsprint( sbuff, tok, s     ) { sbuff+= sprintf( sbuff, " %s: ",s); sbuff+=sprint(sbuff,tok); }
//define _lscan  ( sbuff, tok       ) { sbuff+= sprintf( sbuff, " %s: ",#tok); sbuff+=sprint(sbuff,tok); }
//#define _lsprint_( sbuff, tok, ... ) { _lsprint(sbuff,tok); _lsprint_(sbuff,__VA_OPT__); }
//#define _lsprint_( sbuff, tok, ... ) { _lsprint(sbuff,tok); _lsprint_(sbuff,__VA_ARGS__); } // C-macros are stupid. non recursive, this would not work

#define _toDict( mp, tok ){ mp[#tok] = tok; }
#define _fromDict( mp, tok ){ tok = mp[#tok]; }

inline char * fgetsNonComment(char * str, int num, FILE * stream, char commentChar, int nMaxTry = 100 ){
    char* s = nullptr;
    // no more than 100 comments expected, ensure we don't get stuck in infinite loop
    for (int i = 0; i < nMaxTry; i++) {
        s = fgets(str, num, stream);
        //printf("fgetsNonComment [%i] '%s'", i, s );
        if (s){
            if (s[0] == commentChar)continue;
        }
        break;
    }
    //printf( "fgetsNonComment '%s'", s );
    return s;
}

//  TODO:
// Universal data loading idea:
//  - load all data to std::map<string,string> "craft.velocity"->"0.0 1.0 3.0"
//  - pull named tokens from map and parse them to data   "0.0 1.0 3.0" - > (Vec3d){0.0,1.0,30.}

// TODO: GUI should have also such macros
//       - bind pointers &float to data item to GUI sliders


/*
// making named items in dictionary

#define _toDict  ( char* mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
#define _fromDict( char* mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }

template <typename T> void toBuff  (  const T& v, void* buff, int& i){ (*(T*)(buff+i))=v; i+=sizeof(T); };

#define _toSDict  ( mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
#define _fromSDict( mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
*/


/*
inline char* file2str(char const  *fname ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
    char *buff = NULL;
	fseek (fptr, 0, SEEK_END);
    int length = ftell(fptr);
    buffer = new char[length];
    fseek (fptr, 0, SEEK_SET);
    if (buffer){ fread (buff, 1, length, fptr); }
    fclose (fptr);
	return buff;
}
*/

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
int processFileLines( const char * fname, Func func ){
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

inline int saveBin( const char *fname, int n, char * data ){
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

inline int loadBin( const char *fname, int n, char * data ){
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

inline int fileGetNextKey( FILE  *fptr, const char * keystr, char * tmp ){
	int nch = strlen(keystr);
	//printf( "fileGetNextKey [%s] %i \n", keystr, nch );
    while(fgets(tmp, N_CHAR_TMP, fptr)){
        //printf( "fileGetNextKey: %s %i \n", tmp, strncmp(tmp,keystr,nch) );
        if( strncmp(tmp,keystr,nch)==0 ) return ftell(fptr);
        //printf( "fileGetNextKey: NOOOO! \n" );
    };
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
inline char* fileGetSection(const char *fname, const char * keystr, const char * endstr ){
	//printf( "fileGetSection \n" );
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

inline int whichKey( const char* tmp, int nkey, const char ** keys ){
    for(int ikey=0; ikey<nkey; ikey++){
        const char * key = keys[ikey];
        //printf( "whichKey[%i] %s %s \n", ikey, key, tmp );
        int i=0;
        bool match=true;
        while(key[i]!='\0'){ if(key[i]!=tmp[i])match=false; i++; }
        if(match) return ikey;
    }
    return -1;
}

inline char ** fileGetSections(const char *fname, int nkey, const char ** keys, const char* begstr ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
	char      tmp[N_CHAR_TMP];
	int nb = strlen(begstr);
	char** result = new char*[nkey];
	for(int ikey=0; ikey<nkey; ikey++){ result[ikey] = NULL; }
	int ikey=-1,i0=-1,i1;
	//printf( "fileGetSections\n" );
	while( (i1=fileGetNextKey( fptr, begstr, tmp ))>=0 ){
        //printf( " ikey %i i0  %i i1 %i \n", ikey, i0, i1 );
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
	fclose(fptr);
	return result;
}

inline int checkNullSections( int n, char ** sections ){
    for(int i=0; i<n; i++){ if(sections[i]==NULL) return i; }
    return 0;
}

inline void saveStr( const char * fname, const char * str ){
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
