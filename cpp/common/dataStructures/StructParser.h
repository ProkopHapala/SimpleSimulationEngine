
#ifndef StructParser_h
#define StructParser_h

#include <vector>
#include <cstdio>
#include <cstring>

#include "parsing.h"

namespace parsing{

class StructParser{ public:
    std::vector<ParserItem> items;

    int ich;
    int nch;
    char * str;

    // params
    char cOPEN  = '{';
    char cCLOSE = '}';
    char cNEXT  = ';';

    void parseString( int nch_, char * str_ ){
        str = str_;
        nch = nch_;
        ich = 0;
        //printf( "%s\n", str );
        parse( 0, 0, -1 );
    }

    int parse( int ich0, int level, int parent ){
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
                    items.push_back( (ParserItem){oi,ich-oi,0,0,0, level, parent} );
                    printf("LEAF item[%i] %i %i : %.*s\n", items.size()-1, level, parent,   ich-oi,str+oi );
                }else{
                    items[iCur].nch = ich-oi;
                    //printf( " finished %i \n", ich );
                }
                nbranch++;
                iCur = -1;
                oi   = ich+1;
            //}
            }else if (ch==cOPEN) {
                iCur = items.size();
                items.push_back( (ParserItem){oi,0,ich-oi,0,0,level, parent} );
                printf("NODE item[%i] %i %i : %.*s\n",  items.size()-1,    level, parent,    ich-oi,str+oi );
                parse( ich+1, level+1, iCur );
            }else if (ch==cCLOSE) {
                if( parent>=0 ){
                    items[parent].nbranch = nbranch;
                    items[parent].nitem   = items.size()-parent;
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
        return iit;
    }

    void printItemStruct(){
        for( int i=0; i<items.size(); i++ ){
        //for( ParserItem& item : items ){
            ParserItem& it = items[i];
            //printf("%03i (%04i,%02i) %02i(%02i,%02i) : ", i, it.ibeg, it.iend-it.ibeg, it.level, it.nbranch, it.nitem );
            printf("%03i (%04i,%04i) (%02i,%02i) %03i %03i : ", i, it.ibeg, it.nch,  it.nbranch,it.nitem, it.level, it.parent );
            printf("%*c"   , 4*it.level, ' ' );
            //printf("%.*s\n", it.iend-it.ibeg, str + it.ibeg );
            if( (it.nbranch+it.nitem)>0 ){
                printf("%.*s", it.nch_name, str + it.ibeg );
                printf("|" );
                //printf("%.*s\n", it.iend-it.ibeg-it.nch_name, str + it.ibeg+it.nch_name );
                printf("%.*s\n", it.nch-it.nch_name, str + it.ibeg+it.nch_name );
            }else{
                //printf("%.*s\n", it.iend-it.ibeg, str + it.ibeg );
                printf("%.*s\n", it.nch, str + it.ibeg );
            }
        }
    }

};

}

#endif


