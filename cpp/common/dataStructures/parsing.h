
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

struct ParserNode{
    char kind ;     // in case of bracket this is  "open"-"close"
    char cBond;     // which kind of character bind to the previous token ?
    int  parent;
    int  prev   = -1 ;

    int  ic0;
    int  nc     = -1 ;
    int nbranch = -1 ;
    int level   = -1 ;




    ParserNode() = default;
    ParserNode(char kind, int ic0, int parent, char  ) : kind(kind),ic0(ic0),parent(parent){};

};

}

#endif


