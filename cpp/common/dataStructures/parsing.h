
#ifndef parsing_h
#define parsing_h

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

#endif


