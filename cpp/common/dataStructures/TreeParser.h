
#ifndef TreeParser_h
#define TreeParser_h

#include <vector>
#include <cstdio>
#include <cstring>

#include "parsing.h"


/*

There are 4 kinds of chains with    ( decreasing coupling strength )  <=>  (  increasing split priority )

1)  Name-String                        aa or bracket ()
2)  Void-Binded                        a_a a_() ()_a ()_()     possibly also   a.a a.() ().a ().()
3)  Expressions binded by operators    a+a a+() ()+a ()+()   a+b+c  a*b+c    ...
4)  Tupes with separators              a,b a,() (),a

*/


namespace parsing{

//                                                     !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
char charTypes_1[] = "                                 +`++++`()++;+++aaaaaaaaaa+;+++++aaaaaaaaaaaaaaaaaaaaaaaaaa(+)+aaaaaaaaaaaaaaaaaaaaaaaaaaaa(+)+";
//char charTypes_1[] = "                               -a----a()--,---aaaaaaaaaa-------aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-aaaaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";
//char charTypes[] = "                                 -'----'()--,---0000000000-------aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-a'aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";
//char charPrior_1[] = "                              13 3354 2256 63500000000009 78783000000000000000000000000002 240 000000000000000000000000002423";
char charCoule_1[] = "                                 5 5555   551595          2055555                           5 5                             5 5";
char charCompl_1[] = "                                        )(/- + *            > <                            ] [                             } { ";

struct ParserGlobalState{
/*
    int    level    = 0    ;
    char   cQuote   = 0    ;    // are we in some litaral?      '"`    ?
    char   cBracket = ' '  ;
    char   cPrev    = '('  ;    // previous kind
    char   cChain   = ' '  ;    // type of current chain
    char   chPrior  = ' '  ;    // chain priority
    char   oct      = ' '  ;
    int    parNode  = 0    ;    // parent node
    int    curNode  = 0    ;    // current node
    int    ic       = 0    ;    // position in string
 */

    int    level    = 0    ;
    char   cQuote   = 0    ;    // are we in some litaral?      '"`    ?
    char   cBracket = ' '  ;

    //char   cSuper  = '('   ;    // previous kind
    char   cPrev   = '('   ;    // previous kind
    char   pSuper  = '('   ;    //

    char   iSuper  = -1    ;    // previous kind
    char   iPrev   = -1    ;    // previous kind

    char   oct      = ' '  ;
    int    ic       = 0    ;    // position in string

};







/*
struct ParserLocalState{ // this should be on stack
    int    parent;
    //int    oi;           // previous position for calling recursion

    ParserLocalState(int parent):parent(parent){};
}
*/

class TreeParser{ public:
    std::vector<ParserNode> items;

    int nch;
    char * str;

    ParserGlobalState gs;

    bool bBanEmptyNode       = false;
    bool bEmptyChainOperator = true;


    // params
    /*
    char cOPEN   = '(';
    char cCLOSE  = ')';
    char cNEXT   = ',';
    char cOP     = '-';
    char cQUOTE  = '`';
    */

    //temp
    char ch;
    char ct;

    inline int fail(){
        gs.ic = nch+1;
        return -1;
    }

    inline void reconectTree( int iUp, int iDown ){
        int ipar = items[iDown].parent;
        items[iDown].parent = iUp;
        items[iUp  ].parent = ipar;
    }

    inline void finishNode( ){
        printf( "finishNode \n" );
        items[gs.curNode].nc = gs.ic - items[gs.curNode].ic0;
        int ipar  = items[gs.curNode].parent;
        if(ipar<0){
            printf( "ERROR: str[%i] exit from top-node \n", gs.ic );
            fail();
            return;
        }
        gs.curNode = ipar;
        gs.level--;
    }

    inline void startNode( char kind ){
        printf( "startNode %c ic %i parent %i \n", kind, gs.ic, gs.curNode );
        //items.emplace_back( kind, gs.ic, gs.curNode  );
        items.push_back( ParserNode( kind, gs.ic, gs.iSuper, gs.iPrev ) );
        //items.push_back( ParserNode(kind, gs.ic, gs.curNode, gs.prev ) );
        items.back().level = gs.level;
        gs.curNode=items.size()-1;
        gs.cPrev = kind;
        gs.level++;
    }

    inline void superNode( char kind ){

        int iDown = gs.curNode;

        printf( "superNode %c iDown %i \n", kind, iDown );
        items[iDown].level++;
        finishNode();
        startNode( kind );
        items.back().ic0 = items[iDown].ic0;
        reconectTree( gs.curNode, iDown );
    }

    //inline namedBracket(){
    //    items[curNode].kind='F';
    //}


    inline void nextInChain( char ct, char ch ){
        char coupling = charCoule_1[ch];
        if( coupling < gs.chCoupling ){   // split chains
            superNode(ct);
        }else{                            // next
            startNode(ct);
        }
    }



    int parse( ){
        while( gs.ic < nch ){
            //char ch = str[ich];
            ch = str[gs.ic];         // current character
            ct = charTypes_1[ ch ];  // type of current character
            cp = charPrior_1[ ch ];  // character priority


            printf( "str[%i] %c -> %c \n", gs.ic, ch, ct );

            // ========== quotation
            if( gs.cQuote ){ // in quotation
                if(gs.cQuote==ch){ gs.cQuote=0; }
                gs.ic++; continue;
            }else{        // not int quotatio
                if(ct=='"'){
                    gs.cQuote=ch;
                    gs.ic++; continue;
                }
            }







            // ========== not quotation

            //if (  ct== ' ' ){       // space
            //    gs.bSpace = true;
            //} else


            if ( ct!=gs.ctPrev ){


            }





            if       ( ct== 'a' ) {       // Name
            //    a after a       check previous space, if space    super-node-'_'   (if allowed)      active plit_&_push_down by "_chain" operator
            //    a after )       super-node-'_'  (if allowed)                                         split_&_push_down by "_chain" operator  (if allowed)
            //    a after (       start sub-node
            //    a after +       start sub-node
            //    a after ;       start sub-node
                if ( gs.cPrev=='a' ){
                    if( gs.oct==' ' ){
                        superNode( '_' );
                        startNode( 'a' );
                    }
                }else if ( gs.cPrev==')' ){
                    printf( " )a a_a   cPrev %c \n", gs.cPrev );
                    superNode( '_' );
                    startNode( 'a' );
                }else{
                    printf( " +a ;a  cPrev %c \n", gs.cPrev );
                    startNode( 'a' );         // +a  ;a
                }
            }else if (  ct== '+' ){
            //    + after a       finish a, super-node-'+'
            //    + after )       finish ), super-node-'+'
            //    + after (       start sub-node  unary-operator
            //    + after +       ignore (chain operators)
            //    + after ;       start sub-node  unary-operator
                if( (gs.cPrev == 'a')||( gs.cPrev ==')' ) ){
                    //finishNode();
                    //nextInList( ct, ch );
                    if( gs.cSuper=='+' ){
                        superNode( '+' );
                    }else{
                        superNode( '+' );
                    }
                }else if( gs.cPrev =='(' ){
                    startNode( '~' );
                }else if( gs.cPrev ==';' ){
                    startNode( '~' );
                }
            }else if (  ct== ';'  ){   // separator
            //    ; after a       finish a, super-node-';'
            //    ; after )       finish ), sub-node-';'
            //    ; after (       if(bEmptyNode)  start-sub-node, finish-sub-node,  start sub-node
            //    ; after +       ERROR: dangling +
            //    ; after ;       if(bEmptyNode)  start-sub-node, finish-sub-node,  start sub-node
                if ( (gs.cPrev=='a')||(gs.cPrev==')') ){
                    //finishNode();
                    //nextInList( ct, ch );
                    if( gs.cSuper=='+' ){
                        hyperNode( ';' );
                    }else {
                        superNode( ';' );
                    }
                //}else if (gs.cPrev==')'){
                //    //finishNode();
                //    nextInList( ct, ch );
                }else if (gs.cPrev=='('){
                    startNode ('0');
                    gs.cPrev=';';
                }else if (gs.cPrev=='+'){
                    printf("ERROR: str[%i] dangling operator before separator * %c \n", gs.ic, ch );
                    return fail();
                }else if( (gs.cPrev==';')  ){   // empty node - may be disabled
                    if( bBanEmptyNode && (gs.cPrev==';') ){
                        printf("ERROR: str[%i] chained separator %c, but empty nodes not allowed\n", gs.ic, ch );
                        return fail();
                    }
                    finishNode();
                    startNode (';');
                }
            }else if ( ct== '(' ){    // open bracket
            //    ( after a       Funtion call;    super-node-'_'                   change a-node to a(-node (namedBracket)
            //    ( after )       split_&_push_down by "_chain" operator  (if allowed)
            //    ( after (       start sub-node (
            //    ( after +       start sub-node (
            //    ( after ;       start sub-node (
                if( (gs.cPrev=='a') || (gs.cPrev==')') ){  // named bracket
                    superNode( '_' );
                    startNode( '(' );
                }else{
                    startNode( '(' );
                }
                gs.cBracket = charCompl_1[ch];
            }else if ( ct== ')' ) {        // close bracket
            //    ) after a       finish a, finish bracket
            //    ) after )       finish )
            //    ) after (       finish empty bracket (if allowed)
            //    ) after +       ERROR - usaturated +
            //    ) after ;       finish empty ;, finish bracket
                if( ch != gs.cBracket ){
                    printf( " ERROR: str[%i] closing wrong bracket %c instead of %c \n",  gs.ic, ch, gs.cBracket );
                    return fail();
                }

                if( (gs.cPrev == 'a')||(gs.cPrev ==')') ){
                    finishNode();
                }else if( gs.cPrev =='(' ){
                    if(bBanEmptyNode){
                        printf( " ERROR: str[%i] empty brackets not allowed \n",  gs.ic );
                        return fail();
                    }
                    finishNode();
                }else if( gs.cPrev ==';' ){
                    finishNode();
                }
                gs.cPrev=')';
            }

            gs.oct=ct;
            gs.ic++;
        }
        return 0;
    }

    void parseString( int nch_, char * str_ ){



        str = str_;
        nch = nch_;

        //items.emplace_back( '(', 0, -1, -1 );
        items.push_back( ParserNode( '(', 0, -1 ) );

        gs.ic      = 0   ;
        gs.oct     = ' ' ;
        gs.level   = 0   ;
        gs.curNode = 0   ;
        gs.cQuote  = 0   ;
        gs.cPrev   = '(' ;


        //printf( "%s\n", str );
        parse( );
    }

    void checkTables(){
        char s_[255];
        for(int i=0 ;i<255;i++){ s_[i]=' ';            };
        for(int i=33;i<127;i++){ s_[i]=i;              }; printf("%s\n", s_+32 );
        for(int i=33;i<127;i++){ s_[i]=charTypes_1[i]; }; printf("%s\n", s_+32 );
        for(int i=33;i<127;i++){ s_[i]=charCompl_1[i]; }; printf("%s\n", s_+32 );
    }

    void printItems(){
        printf("=============\n");
        for( int i=0; i<items.size(); i++ ){
            const ParserNode& it = items[i];
            printf( "it[%i]    %.*s \n",i,     it.nc, str+it.ic0 );
        }
    }



    /*

        inline void nextInList( char ct, char ch ){

    }

    int parse( ){
        while( gs.ic < nch ){
            //char ch = str[ich];
            ch = str[gs.ic];         // current character
            ct = charTypes_1[ ch ];  // type of current character

            printf( "str[%i] %c -> %c \n", gs.ic, ch, ct );

            // ========== quotation
            if( gs.cQuote ){ // in quotation
                if(gs.cQuote==ch){ gs.cQuote=0; }
                gs.ic++; continue;
            }else{        // not int quotatio
                if(ct=='"'){
                    gs.cQuote=ch;
                    gs.ic++; continue;
                }
            }

            // ========== not quotation

            //if (  ct== ' ' ){       // space
            //    gs.bSpace = true;
            //} else

            if       ( ct== 'a' ) {       // Name
            //    a after a       check previous space, if space    super-node-'_'   (if allowed)      active plit_&_push_down by "_chain" operator
            //    a after )       super-node-'_'  (if allowed)                                         split_&_push_down by "_chain" operator  (if allowed)
            //    a after (       start sub-node
            //    a after +       start sub-node
            //    a after ;       start sub-node
                if ( gs.cPrev=='a' ){
                    if( gs.oct==' ' ){
                        superNode( '_' );
                        startNode( 'a' );
                    }
                }else if ( gs.cPrev==')' ){
                    printf( " )a a_a   cPrev %c \n", gs.cPrev );
                    superNode( '_' );
                    startNode( 'a' );
                }else{
                    printf( " +a ;a  cPrev %c \n", gs.cPrev );
                    startNode( 'a' );         // +a  ;a
                }
            }else if (  ct== '+' ){
            //    + after a       finish a, super-node-'+'
            //    + after )       finish ), super-node-'+'
            //    + after (       start sub-node  unary-operator
            //    + after +       ignore (chain operators)
            //    + after ;       start sub-node  unary-operator
                if( gs.cPrev == 'a' ){
                    //finishNode();
                    superNode( '+' );
                }else if( gs.cPrev ==')' ){
                    //finishNode();
                    superNode( '+' );
                }else if( gs.cPrev =='(' ){
                    startNode( '~' );
                }else if( gs.cPrev ==';' ){
                    superNode( '~' );
                }
            }else if (  ct== ';'  ){   // separator
            //    ; after a       finish a, super-node-';'
            //    ; after )       finish ), sub-node-';'
            //    ; after (       if(bEmptyNode)  start-sub-node, finish-sub-node,  start sub-node
            //    ; after +       ERROR: dangling +
            //    ; after ;       if(bEmptyNode)  start-sub-node, finish-sub-node,  start sub-node
                if (gs.cPrev=='a'){
                    //finishNode();
                    superNode(';');
                }else if (gs.cPrev==')'){
                    finishNode();
                    startNode(';');
                }else if (gs.cPrev=='('){
                    startNode ('0');
                    gs.cPrev=';';
                }else if (gs.cPrev=='+'){
                    printf("ERROR: str[%i] dangling operator before separator * %c \n", gs.ic, ch );
                    return fail();
                }else if( (gs.cPrev==';')  ){   // empty node - may be disabled
                    if( bBanEmptyNode && (gs.cPrev==';') ){
                        printf("ERROR: str[%i] chained separator %c, but empty nodes not allowed\n", gs.ic, ch );
                        return fail();
                    }
                    finishNode();
                    startNode (';');
                }
            }else if ( ct== '(' ){    // open bracket
            //    ( after a       Funtion call;    super-node-'_'                   change a-node to a(-node (namedBracket)
            //    ( after )       split_&_push_down by "_chain" operator  (if allowed)
            //    ( after (       start sub-node (
            //    ( after +       start sub-node (
            //    ( after ;       start sub-node (
                if( (gs.cPrev=='a') || (gs.cPrev==')') ){  // named bracket
                    superNode( '_' );
                    startNode( '(' );
                }else{
                    startNode( '(' );
                }
                gs.cBracket = charCompl_1[ch];
            }else if ( ct== ')' ) {        // close bracket
            //    ) after a       finish a, finish bracket
            //    ) after )       finish )
            //    ) after (       finish empty bracket (if allowed)
            //    ) after +       ERROR - usaturated +
            //    ) after ;       finish empty ;, finish bracket
                if( ch != gs.cBracket ){
                    printf( " ERROR: str[%i] closing wrong bracket %c instead of %c \n",  gs.ic, ch, gs.cBracket );
                    return fail();
                }

                if( (gs.cPrev == 'a')||(gs.cPrev ==')') ){
                    finishNode();
                }else if( gs.cPrev =='(' ){
                    if(bBanEmptyNode){
                        printf( " ERROR: str[%i] empty brackets not allowed \n",  gs.ic );
                        return fail();
                    }
                    finishNode();
                }else if( gs.cPrev ==';' ){
                    finishNode();
                }
                gs.cPrev=')';
            }

            gs.oct=ct;
            gs.ic++;
        }
        return 0;
    }


    */





    /*

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

    */

};


}; // namespace parsing

#endif


