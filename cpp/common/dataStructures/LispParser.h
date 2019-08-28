
#ifndef LispParser_h
#define LispParser_h

#include <vector>
#include <cstdio>
#include <cstring>

#include "parsing.h"

// inspired by http://zserge.com/jsmn.html and https://github.com/zserge/jsmn

class LispParser{ public:
    std::vector<ParserItem> items;

    int ich;
    int nch;
    char * str;

    // params
    char cOPEN  = '(';
    char cCLOSE = ')';
    char cNEXT  = ',';

    void parseString( int nch_, char * str_ ){
        str = str_;
        nch = nch_;
        ich = 0;
        //printf( "%s\n", str );
        parse( 0, 0, nullptr );
    }

    int parse( int ich0, int level, ParserItem* parent ){
        int oi      = ich0;
        ich         = oi;
        int nbranch =  0;
        int iCur    = -1;
        //printf( "%s\n", str );
        //printf( "ich0 %i level %i \n", ich0, level );
        while( ich < nch ){
            char ch = str[ich];
            //printf( "%i : %c\n", ich, ch );
            if (ch==cNEXT){
                if( iCur<0 ){
                    items.push_back( (ParserItem){oi,ich,0,0,0, level} );
                    //printf(" item[%i] (%i,%i) : %.*s\n", items.size(), oi, ich, ich-oi, str+oi );
                }else{
                    items[iCur].iend = ich;
                    //printf( " finished %i \n", ich );
                }
                nbranch++;
                iCur = -1;
                oi   = ich+1;
            //}
            }else if (ch==cOPEN) {
                iCur = items.size();
                items.push_back( (ParserItem){oi,0,ich-oi,0,0} );
                //printf(" item[%i] (%i,%i) : %.*s\n", items.size(), oi, ich, ich-oi, str+oi );
                parse( ich+1, level+1, &items[iCur] );
            }else if (ch==cCLOSE) {
                if( parent ){
                    parent->nbranch = nbranch;
                    parent->nitem   = items.size()-( parent - &items[0] );
                }
                return 0;
            }
            ich++;
        }
        return 0;
    }

    int get( const char * name ){
        // use like:  parser.get( "Name1.Name3.Name6" );
        // fast search for named item within string
        const char *  id  = name;
        int nch     = strchr( id, '.') - id;
        int iit = 0;
        while(iit < items.size() ){
            ParserItem& it = items[iit];
            if( strncmp( id, str+it.ibeg,  nch ) == 0 ){ // found
                id  += nch;
                nch  = strchr( id, '.') - id;
            }else{
                iit += it.nitem;  // skip all sub-items
            }
        }

    }

    void printItemStruct(){
        for( int i=0; i<items.size(); i++ ){
        //for( ParserItem& item : items ){
            ParserItem& it = items[i];
            printf("%06i (%i,%i) : ", i, it.ibeg, it.iend );  printf("%*c", 4*it.level, ' ' );  printf("%.*s\n", it.iend-it.ibeg, str + it.ibeg );
        }
    }

};

#endif


