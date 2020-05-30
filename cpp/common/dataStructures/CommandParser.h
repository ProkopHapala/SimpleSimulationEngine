
#ifndef CommandParser_h
#define CommandParser_h

#include <stdlib.h>

// simple command parser
// # comment
// commands are separated by some separator (;, \t )
// first world is name of command
// each application overloads its own CallTable function

//#include <vector>
#include <unordered_map>
#include <cstdio>
#include <cstring>

// TODO : dynamic construction of variables is quite complicated, for now we don't do it
//        this means that commands can have just constant arguments of primitive type
class VarObj{
    char  type;
    void * prt = NULL;
};

class FunctionObj{
    //char name    []
    char argTypes[16];

    //void exec( char* line,  int nToks, int * tokIs, unordered_map<string,VarObj> varDict ){
    void exec( char* line,  int nToks, int * tokIs ){
        for( int i=0; i<(nToks-1); i++ ){
            char type = argTypes[i];
            char * tok = line + tokIs[i+1];
            switch(type){
                case  'f': break; // float
                case  'i': break; // int
                case  's': break; // string

            }
        }

    }
};

//================== test commands

void command1( int   a, float b        ){ printf(" command1 ( %i, %f ) \n", a, b );     }
void command2( float a, float b, int c ){ printf(" command2 ( %f, %f, %i ) \n", a, b, c ); }
void command3( float a, char* s, int c ){ printf(" command3 ( %f, %s, %i ) \n", a, s, c ); }

//==================

class CommandParser{ public:
    static const int maxLineLen = 2048;
    static const int maxTok     = 128;
    int       maxLines   = 10000;
    char comentChar='#';
    char varChar   ='$';
    //vector<int> tokIs;
    int nToks = 0;
    //int tokIs[maxTok];
    char* toks[maxTok];
    char  line[maxLineLen];

    int iline,iarg;

    void line2toks( char* line ){
        //tokIs.clear();
        for(int i=0; i<nToks; i++) toks[i]=NULL; // clean previous tokens
        nToks=0;
        bool inString = false;
        bool inTok    = false;
        for(int i=0; i<maxLineLen; i++ ){
            if  (  (!inString) &&  ((line[i]==' ' )||(line[i]=='\t'))){
                line[i]='\0';
                inTok=false;
            }else if ( (line[i]=='\0')||(line[i]==comentChar)){
                break;
            }else if (!inTok){
                if (line[i]=='"'){ inString = true; i++; }
                //tokIs.pushBack(i);
                inTok=true;
                //tokIs[nToks]=i;
                toks[nToks]=line+i;
                nToks++;
            }else{
                if (line[i]=='"'){ inString = false; line[i]='\0'; }
            };
        }
    }

    int execFile( char * fname ){
        FILE *pFile=fopen( fname,"r");
        if(pFile==NULL){
            printf("Failed to load %s \n", fname );
            return -1;
        }
        for(iline=0; iline<maxLines; iline++){
            char *str = fgets( line, maxLineLen, pFile );
            if( str==NULL ) break;
            if( str[0]==comentChar ) continue;
            //printf( "%i : %s", iline, str );
            //exit(0);
            line2toks( str );
            iarg = 0;
            callTable( toks[0] );
        };
        return 0;
    }

    void argError(){ printf( "ERROR line %i arg %i", iline+1, iarg );  exit(0); };

    int   _int  (char* s){ int   i; iarg++; if(1!=sscanf( s, "%i", &i )) argError(); return i; };
    float _float(char* s){ float f; iarg++; if(1!=sscanf( s, "%f", &f )) argError(); return f; };

    void callTable_default( char * funcName ){
        //char * funcName = line + tokIs[0];
        //printf( "funcName : %s \n", funcName );
        // TODO: strcmp should be replaced by hashtable
        if      ( 0==strcmp(funcName, "command1" ) ){ command1( _int  (toks[1]), _float(toks[2])                ); }
        else if ( 0==strcmp(funcName, "command2" ) ){ command2( _float(toks[1]), _float(toks[2]), _int(toks[3]) ); }
        else if ( 0==strcmp(funcName, "command3" ) ){ command3( _float(toks[1]),        toks[2] , _int(toks[3]) ); };
    }

    virtual void callTable( char * funcName ){
        callTable_default( funcName );
    }

};

#endif


