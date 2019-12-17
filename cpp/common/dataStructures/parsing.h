
#ifndef parsing_h
#define parsing_h

namespace parsing{
//
//static const char ascii_mask_[] = "";


/*
void bakeAsciiMask(const char* cases, bool* mask){
    for(int i=0;i<256;i++){
        char* pc=cases;
        while((*pc)>0){
            if((*pc)==i){ mask[i]=true; }else{  mask[i]=false; }
            pc++;
        }
    }
}
*/

//                                                   !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
char charTypes[] = "                                 -'----'()--,---0000000000-------aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-a'aaaaaaaaaaaaaaaaaaaaaaaaaa(-)-";

struct ParserItem{
    int ibeg;
    int nch;      // number of charactrs
    int nch_name; // number of characters before '{'
    int nbranch;  // number of direct branches
    int nitem;    // total number of other items contained
    int level;
    int parent;
    //int type;
};

struct Node{
    char kind ;     // in case of bracket this is  "open"-"close"
    char cBond;     // which kind of character bind to the previous token ?
    int  parent;
    int  prev   = -1 ;

    int  ic0;
    int  nc     = -1 ;
    int nbranch = -1 ;
    int level   = -1 ;

    Node() = default;
    Node(char kind, int ic0, int parent, char  ) : kind(kind),parent(parent),ic0(ic0){};

};


struct Tuple{
    char bra   ='_';   //   ( { [ _
    char cKind;  //   a o _ + ;
    //char pSep;   //   separator priority

    int  iPrev  = -1 ;
    int  iSuper = -1 ;

    int  ic;
    int  nc     = -1 ;
    int nbranch = -1 ;
    int level   = -1 ;

    Tuple() = default;

    Tuple(char cKind, char bra, int ic, int iSuper ) : bra(bra),cKind(cKind),iSuper(iSuper),ic(ic){
        //bra     = '_';
        iPrev   = -1 ;
        iSuper  = -1 ;
        nc      = -1 ;
        nbranch = -1 ;
        level   = -1 ;
    };

    void print(){ printf( "Tuple{bra`%c` kind`%c` iSup%i lev%i ic%i}",   bra,  cKind,  iSuper, level, ic); };

    //Tuple(char cKind, char cBracket, int ic, int iSuper, int iPrev  ) : cKind(cKind),bra(bra),ic(ic),iSuper(iSuper),iPrev(iPrev){};

};

struct Break{
    int  ic;
    char cl,cr;
    char cb;
    Break(int ic, char cl, char cr, char cb ): ic(ic), cl(cl), cr(cr), cb(cb) {}
};

}

#endif


