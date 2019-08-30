
#ifndef TreeParser_h
#define TreeParser_h

#include <vector>
#include <cstdio>
#include <cstring>

#include "macroUtils.h"
#include "parsing.h"


/*

There are 4 kinds of chains with    ( decreasing coupling strength )  <=>  (  increasing split priority )

1)  Name-String                        aa or bracket ()
2)  Void-Binded                        a_a a_() ()_a ()_()     possibly also   a.a a.() ().a ().()
3)  Expressions binded by operators    a+a a+() ()+a ()+()   a+b+c  a*b+c    ...
4)  Tupes with separators              a,b a,() (),a

*/


namespace parsing{


enum class TupleKind {
name=0,
op  =1,
expr=3,
lst =9,
};




//                                                      !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
char charTypes_1[]  = "                                 +`++++`()++;+++aaaaaaaaaa+;+++++aaaaaaaaaaaaaaaaaaaaaaaaaa(+)+aaaaaaaaaaaaaaaaaaaaaaaaaaaa(+)+";
//char charTypes_1[]  = "                               -a----a()--,---aaaaaaaaaa-------aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-aaaaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";
//char charTypes[]  = "                                 -'----'()--,---0000000000-------aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-a'aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";
//char charPrior_1[]  = "                                13 3354 2256 63500000000009 78783000000000000000000000000002 240 000000000000000000000000002423";
char charPrior_1[]  = "                                01 1111   1181110000000000191111100000000000000000000000000   13 00000000000000000000000000 1 1";
//char charCouple_1[] = "                                 5 5555   551595          2055555                           5 5                             5 5";
char charCompl_1[]  = "                                        )(/- + *            > <                            ] [                             } { ";

struct ParserGlobalState{
    //Tuple  tup;
    char   cQuote   = 0     ;  // are we in some litaral?      '"`    ?
    //char   cBracket = ' ' ;
    //char   cMode    = ' ' ;
    char   cPrev    = ' '   ;  // _ a +
    char   ctSep    = ' '   ;  // _ ; ( )    previous separator
    char   chSep    = ' '   ;
    //int   iSuper  = -1    ;   // previous kind
    //int   iPrev   = -1    ;   // previous kind
    char    icur    = 0    ;    // current tuple
    char   oct      = ' '  ;
    int    ic       = 0    ;    // position in string

    //int level  = 0;
};

class TreeParser{ public:
    //std::vector<ParserNode> items;
    //std::vector<Break> breaks;

    std::vector<Tuple> items;

    int nch;
    char * str;

    ParserGlobalState gs;

    //bool bBanEmptyNode       = false;
    bool bBanEmptyNode       = true;
    bool bEmptyChainOperator = true;
    bool bBanRightDanglingOperator = true;

    //temp
    char ch;
    char ct;
    char cp;

    inline int fail(){
        gs.ic = nch+1;
        return -1;
    }

    inline void changeTuple( char cKind ){ items[gs.icur].cKind = cKind; }

    inline void startTupleBra( char cKind, char bra  ){
        //printf( "startTuple  (%c |%c|%i)", cKind, bra,bra==' ' );
        //int level = 0; if(gs.icur>0){level=items[gs.icur].level + 1;}

        items.push_back( Tuple( cKind, bra, gs.ic, gs.icur ) );
        //items.back().level = gs.level;
        gs.icur = items.size()-1;
        printf( "new[%i] ", gs.icur ); items.back().print(); puts(" ");
        //gs.level++;

    }

    inline void startTuple( char cKind  ){ startTupleBra( cKind, '_'  ); }

    bool endTuple( ){
        //printf( "finishNode \n" );
        Tuple& it = items[gs.icur];
        if( bBanEmptyNode && (it.cKind==' ') ){ printf( "ERROR str[%i]: tuple [%i] ended empty \n", gs.ic ); exit(0); return false; }
        it.nc = gs.ic - it.ic;
        int iSup  = it.iSuper;
        if(iSup<0)                            { printf( "ERROR str[%i]: exit from top-node \n",     gs.ic ); exit(0); return false; }
        //printf( "endTuple (( it[%i] cK %c   `%.*s` ))   %i->%i \n",gs.icur,  it.cKind,  it.nc, str+it.ic,  gs.icur, iSup );
        printf( "end[%i] ", gs.icur ); items.back().print(); puts(" ");
        gs.icur = iSup;
        //gs.level--;
        return true;
    }

    inline void reconectWrap( int iUp, int iDown ){
        printf( "reconect %i under %i \n", iDown, iUp );
        int iSup = items[iDown].iSuper;
        printf( "Dw[%i].Sup=%i   Up[%i].Sup=%i \n",  items[iDown].iSuper,iUp,   iUp, iSup );
        items[iUp  ].iSuper = iSup;
        items[iDown].iSuper = iUp;

        items[iUp  ].ic = items[iDown].ic;

        //int lev = items[iDown].level;
        //items[iUp  ].level = lev;
        //items[iDown].level = lev+1;

        _swap( items[iUp  ].bra, items[iDown ].bra );

        printf( "iUp[%i] ", iUp   ); items[iUp  ].print(); puts(" ");
        printf( "iDw[%i] ", iDown ); items[iDown].print(); puts(" ");
    }

    //void wrapTuple( char cKind ){
    //    //printf( "(wrapTuple %c)", cKind );
    //    int iDown = gs.icur;
    //    endTuple();
    //    startTuple( cKind );
    //    reconectWrap( gs.icur, iDown );
    //}

    //void wrapTupleAndSub( char cKind, char ct ){
    //    wrapTuple ( cKind   );
    //    startTuple( ct );
    //}



    // ---- DOES NOT WORK - because child nodes are not reconnected
    //void pseudoWrap( int iDown, int cWant ){
    //    //startTuple( cWant );
    //    //reconectWrap( gs.icur, iDown );
    //    char ocK = items[iDown].cKind;
    //    items[iDown].cKind = cWant;
    //    //startTuple( ocK );
    //    int ic = items[iDown].ic;
    //   gs.icur = iDown;
    //    items.push_back( Tuple( ocK, '_', ic, gs.icur ) );
    //    items.back().nc =  gs.ic - ic;
    //}

    void finishLastRec(){
        while( items[gs.icur].level != 0 ){
            endTuple();
        }
    }

    void treeInsert( char cWant ){
        char pwant  = charPrior_1[ cWant ];
        char cKind  = items[gs.icur].cKind;
        char phave  = charPrior_1[ cKind ];
        int  iDown  = -1;
        int levelUp=0;
        while( pwant > phave ){
            iDown  = gs.icur;
            if( items[gs.icur].bra != '_' ){
                printf( "item[%i]HIT-BRAKET|%c|%i(pwant %i(%c)have %i(%c))\n", gs.icur, items[gs.icur].bra, items[gs.icur].bra=='_', pwant-'0',cWant, phave-'0',cKind );
                //endTuple();
                startTuple( cWant );
                reconectWrap( gs.icur, iDown );
                //pseudoWrap( iDown, cWant );
                return;
            }
            //printf( "while[%i] icur %i bra |%c| %i \n", levelUp, gs.icur, items[gs.icur].bra, ' '==items[gs.icur].bra ); if(levelUp>5) exit(0);
            endTuple();
            levelUp++;
            cKind  = items[gs.icur].cKind;
            phave  = charPrior_1[ cKind ];
        }
        printf( "treeInsert lUp %i pwant %i(%c) phave %i(%c) \n", levelUp, pwant-'0',cWant, phave-'0',cKind );
        if( pwant < phave  ){ // we got higher than we wanted => make sub node
            startTuple( cWant );
            if(iDown>0) reconectWrap( gs.icur, iDown );
            //if(iDown>0){ pseudoWrap( iDown, cWant ); }
            //else       { startTuple( cWant );        }
        }
    }

    void putUnderCommon( char cSuper, char ct, char bra ){
        treeInsert( cSuper );
        startTupleBra( ct, bra );
    }

    //void endTupleRecur(char ch, char ct){
    //    char prior  = charPrior_1[ ch ];
    //    char cKind  = items[gs.icur].cKind;
    //    char cPrior = charPrior_1[ cKind ];
    //    int  iDown  = -1;
    //    while( prior > cPrior ){
    //       //wrapTuple();
    //        iDown  = gs.icur;
    //        endTuple();
    //        cKind  = items[gs.icur].cKind;
    //        cPrior = charPrior_1[ cKind ];
    //    }
    //    // we know that prior <= cPrior
    //
    //    //if( prior == cPrior ){   //   (  ,  (()+()) ,  )
    //    //    addSub(' ');
    //    //}else{
    //    //    wrapTuple(iclosed);  //      ;( (()+()) ,  )
    //    //}
    //    startTuple( ct );
    //    if(iDown>0) reconectWrap( gs.icur, iDown );
    //}

    bool closeBracket(char ket){
        char bra  = charCompl_1[ket];
        char cBra = items[gs.icur].bra;
        while ( cBra=='_' ){
            endTuple();
            cBra = items[gs.icur].bra;
        }
        if( bra == cBra ){
            endTuple();
            return true;
        }
        return false;
    }


    int parse(){
        for(gs.ic=0; gs.ic<nch; gs.ic++ ){

            ch = str        [ gs.ic ]; // current character
            ct = charTypes_1[ ch    ]; // type of current character

            //  -   for the moment   we create formula-kind of tuples with alternating   aaaa_++++_BBB   parts
            char cKind = items[gs.icur].cKind;
            char ctK   = charTypes_1[ cKind ];

            printf( "====== parse[%i]  %c : %c    cK %c  iCur %i  \n",   gs.ic,  ch, ct,  cKind,  gs.icur     );

            if       ( ct=='a'){

                if ( (cKind == '_') || (ctK == ';') ){
                    startTuple( 'a' );
                }else if ( (gs.cPrev ==')') || (cKind=='+') || ( (cKind=='a') && (gs.oct==' ') ) ){
                    //wrapTupleAndSub( '_', ct );       printf(" wrapAndSub   ");
                    putUnderCommon( '_', ct, '_' );
                }else if ( cKind==' ' ){
                    changeTuple( 'a' ); printf( "changeTuple\n" );
                }

            }else if( ct=='+'){

                if ( (cKind == '_') || (ctK == ';') ){
                    startTuple( '+' );
                }else if (  (gs.cPrev ==')') || (ctK=='a' ) ){
                    //wrapTupleAndSub( '_', ct );       printf(" wrapAndSub   ");
                    putUnderCommon( '_', ct, '_' );
                }else if ( cKind==' ' ){
                    changeTuple( 'a' );
                }

            }else if( ct==';' ){

                if     ( cKind == '+' ){ if(bBanRightDanglingOperator){ printf( "ERROR str[%i]: operator before   spearator %c \n", gs.ic, ch ); return -1;  }}
                else if( cKind == ' ' ){ if(bBanEmptyNode            ){ printf( "ERROR str[%i]: empty node before spearator %c \n", gs.ic, ch ); return -1;  }}
                //endTupleRecur(ch, ct );
                treeInsert( ch );
            }else if( ct=='(' ){

                if ( (cKind == '_') || (cKind == ';') ){
                    startTupleBra( ' ', ch );
                } else if( (gs.cPrev ==')') || (cKind=='+') || (cKind=='a' )){
                    //wrapTupleAndSub( '_', ' ' );
                    printf( "ct=='('  to putUnderCommon   \n" );
                    putUnderCommon( '_', ' ', ch );
                    //treeInsert( ch );
                    items[gs.icur].bra = ch;
                }

            }else if( ct==')'){
                if( gs.cPrev == '+'   ){ printf( "ERROR str[%i]: operator before bracket %c \n",             gs.ic, ch                     ); return -1; }
                if( !closeBracket(ch) ){ printf( "ERROR str[%i]: closing bracket %c while %c is opended \n", gs.ic, ch, items[gs.icur].bra ); return -1; }
            }

            printf( "\n "   );

            gs.oct = ct;
        }
        finishLastRec();
        return 0;
    }


    void parseString( int nch_, char * str_ ){
        str   = str_;
        nch   = nch_;
        //printf( "%s\n", str );

        //items.reserve( 10 + (nch/20) );



        //Tuple  tup;
        //gs.cQuote   = 0     ;  // are we in some litaral?      '"`    ?
        gs.cPrev    = ' '   ;  // _ a +
        gs.icur     = -1    ;    // current tuple
        gs.oct      = ' '   ;

        gs.ic       = 0     ;    // position in string
        //gs.level    = 0     ;

        startTupleBra( ';', '(' );
        items.back().level = 0;

        parse();
    }

    int evalLevel(int i){
        int level = 1;
        int j = items[i].iSuper;
        //printf("===== %i \n", i );
        //printf( "[%i]  %i->%i L=%i | %i \n", i, i, j, level, items[j].level );
        while( items[j].level < 0 ){
            //printf( "[%i]  %i->%i L=%i | %i \n", i, j, items[j].iSuper, level, items[j].level );
            j = items[j].iSuper;
            level++;
        }
        //printf("level %i \n", level);
        level += items[j].level;
        items[i].level = level;
        j = items[i].iSuper;
        for(int i=0;i<level;i++){
            level--;
            items[j].level = level;
            j = items[j].iSuper;
        }
    }

    void evalAllLevels(){
        for(int i=1;i<items.size();i++){  evalLevel(i);  }
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
            const Tuple& it = items[i];
            //printf("%*c"   , 4*it.level, ' ' );
            //printf( "it[%i]@%i `%c`%i %c    `%.*s` \n",i,it.iSuper,    it.bra,it.bra=='_',  it.cKind,  it.nc, str+it.ic );
            printf( "it[%02i]@%i `%c`%i %c    ",i,it.iSuper,    it.bra,it.bra=='_',  it.cKind );
            //printf( "%*c", 4*it.level, ' ' );

            printf( "L[%02i]%*c  ", it.level, 4*it.level );
            //printf("\033[0;36m");
            //printf("\033[0;%02im",  i+30   );
            printf( "`%.*s` \n",  it.nc, str+it.ic );
            printf("\033[0m");
        }
    }



    void printLevelColor(){
        int nMaxLevel=10;
        for(int ilev=0; ilev<nMaxLevel; ilev++){
            int ic = 0;
            for(int i=0;i<items.size();i++){
                const Tuple& it = items[i];
                if( it.level==ilev){
                    int ncleft = (it.ic+1)-ic;
                    if(ncleft>0){
                        printf("%*c", ncleft );
                        ic+=ncleft;
                    }
                    printf("\033[0;%02im",  31+(i%8) );
                    printf( "%.*s",  it.nc, str+it.ic );
                    ic+=it.nc;
                }
            }
            printf("\033[0m");
            printf("\n");
        }
    }





};


}; // namespace parsing

#endif


