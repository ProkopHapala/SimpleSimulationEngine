
#ifndef AST_h
#define AST_h

#include "macroUtils.h"
#include "parsing.h"

namespace Lang{

using namespace parsing;

struct Node{
    //int parsing;
    ParserItem* token=0;
    //Node* parent =  0;
    int iparent  = -1;
    int isub0 = -1;
    int nsub  =  0;
    int level =  0;
    int ith   = -1;  // number in current scope
};

class AST{ public:
    int n;
    ParserItem* tokens;
    Node * nodes = 0;
    int  * n2sub = 0;

void reallocOwn(int n){
    // "items" is typically pointer
    _realloc(nodes,n);
    _realloc(n2sub,n);
};

void items2ast(){
    Node* nodes = new Node[n];
    for(int i=0; i<n; i++){
        ParserItem& it = tokens[i];
        Node&       nd = nodes[i];
        int ipar  = tokens[i].parent;
        nd.token  = &it;
        nd.iparent = ipar;
        nd.level  = it.level;
        nodes[ipar].nsub++;
    }
}

void makeParent2subMap(){
    int isub=0;
    for(int i=0; i<n; i++){ // second pass to count isub0
        int& ns=nodes[i].nsub;
        if(ns<=0) continue;
        nodes[i].isub0=isub;
        isub+=ns;
        ns=0;
    }
    for(int i=0; i<n; i++){
        int ipar  = tokens[i].parent;
        int& ns   = nodes[ipar].nsub;
        int  im   = nodes[ipar].isub0 + ns;
        n2sub[im] = i;
        ns++;
    }
}

Node* getSub(int inode, int isub){
    if( nodes[inode].nsub>isub ){
        return nodes+n2sub[ nodes[inode].isub0 + isub ];
    }
    return 0;
}


};




};

#endif


