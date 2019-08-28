
#ifndef AST2Lang_h
#define AST2Lang_h

#include "parsing.h"
#include "AST.h"

#include <unordered_map>

namespace Lang{


struct Symbol{
    std::string name;
};

struct Type : Symbol{
    std::string name;
};

struct Variable : Symbol{

};

struct Function : Symbol{
    std::string name;
};







struct LocalName {
    //std::string name;
    int ith; // from which item in the scope it was defined ?
    //int kind;
    int id;
};

struct LocalNames{
    // for each node we define local namespace, i.e. names which were declared in this scope
    std::unordered_map<std::string,LocalName> names;
};


class AST2L{
    char* code_str = 0;
    AST* ast                    = 0;
    LocalNames**         locals = 0;
    std::vector<Symbol>  symbols;
    //std::vector<>  locals;

int getKnownSymbol(int inode, const std::string& symbol){
    Node* nodes = ast->nodes;
    Node& nd    = nodes[inode];
    int ipar    = nd.iparent;

    while( ipar>=0 ){

        LocalNames* loc = locals[ipar];
        if( loc ){
            auto got = loc->names.find(symbol);
            if(got == loc->names.end() ) continue;
            if( got->second.ith < nd.ith ){
                return got->second.id;
            }
        }
        ipar = nodes[ipar].iparent;
    }
    return -1;
}


void insertSymbols(){
    for(int i=0;i<ast->n;i++){
        ParserItem& tok = ast->tokens[i];

        char* sname = code_str+tok.ibeg;
        std::string name;
        int id       = getKnownSymbol(i, name);
        bool bKnown  = (id>=0);
        bool bCreate = sname[0]!='$';
        if(bCreate){
            //if(bKnown)
            //int id = ;
            int id = symbols.size();
            symbols.push_back( Symbol{ name } );
            //symbols.back().
            if( !locals[i] ) locals[i]=new LocalNames();
            locals[i]->names.insert( {name, (LocalName){ast->nodes[i].ith,id} }  );
        }else{

        }
    }
}

};


/*
std::unordered_map<std::string,int> mKeyws;

std::unordered_map<std::string,int> mVars;
std::unordered_map<std::string,int> mFuncs;
std::unordered_map<std::string,int> mTypes;

std::vector<Keyword>  keyws;
std::vector<Variable> vars;
std::vector<Function> funcs;
std::vector<Type>     types;
*/



/*



class sparseArray{ };


template<class T>
struct MapedArray{
    std::unordered_map<std::string,int> dict;
    std::vector       <T>               arr;

    T* getByName(std::string s){
        auto got = mymap.find(input);
        if ( got == mymap.end() ){ return 0;                 }
        else                     { return &arr[got->second]; }
    }

}

MapedArray keyws;
MapedArray types;
MapedArray funcs;
MapedArray vars;
*/

/*
struct Label{
    std::string s;
    std::unordered_map<std::string,int> ;
};
*/


}


#endif


