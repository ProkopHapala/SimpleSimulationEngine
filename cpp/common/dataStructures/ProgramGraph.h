
#ifndef ProgramGraph_h
#define ProgramGraph_h

#include "parsing.h"
#include "AST.h"
#include <vector>
#include <unordered_map>
#include <string>

/*

Hasing scoped names
 - each named term (Type,Function,Variable) is defined in particular scope. This scope has particular hash. Term is identified not only by it's name, but also by its hash

Variable definition:
    $x:Integer 10;
Function definition:
    $!sqrt <x:Real >y:Real { ;;; };
Type definition:
    $@Complex:Number { x:Real; y:Real; };

Variable assignment
x 10;
Anonymou function
!(<x,>y){y=x*x;}


 for-cyctle is function which broadcast other function
 for( $i:int, 10, !{y=i*i;} );
 for( i, !(>y){y=i*i;} );



*/

unsigned long str_hash( int n, unsigned char* str, unsigned long  h = 5381){
    // from https://stackoverflow.com/questions/7666509/hash-function-for-string
    // http://www.cse.yorku.ca/~oz/hash.html
    //unsigned long  h    = 5381;
    unsigned char* last = str+n;
    for(unsigned char* c=str;c<last;c++) h = ((h << 5) + h) + *c;
    //while (c = *str++) h = ((h << 5) + h) + c; // hash * 33 + c
    return h;
}


namespace ProgramGraph{

enum class TokKind{ none,type,func,var,arg,literal };

//typedef uchar unsigned char;
//uint8_t


char ch_opSep  =';';
char ch_tokSep =' ';
char ch_def    ='$';
char ch_Func   ='!';
char ch_Type   ='#';

unsigned char* source_code = 0; /// whole source code - you can search names inside

struct Scope;
struct CodeBlock;

struct Term{
    Scope* scope=0;            /// in which scope it is defined ?
    unsigned long   nameHash; /// hash of the name for fast search
    //int    id;     /// unique identificator
    int    isrc;     /// index of name in source code
    int    nsrc;     /// number of characters of name in source code

    long makeHash(){
        //nameHash = str_hash( nsrc, source_code+isrc, scope.nameHash );
        nameHash = str_hash( nsrc, source_code+isrc );
        return nameHash;
    }

    /*
    bool equals(const Term& t){
        if(nsrc!=t.nsrc) return false;
        if(scope!=t.scope) return false;
        for(int i=0;i++ nsrc){
            if(source_code[isrc+i]!=source_code[t.isrc+i]) return false;
        }
        return true;
    }
    */

    bool equalName(int i0, int n)const{
        if(nsrc!=n) return false;
        //if(scope!=t.scope) return false;
        for(int i=0;i<nsrc;i++){
            if(source_code[isrc+i]!=source_code[i0+i]) return false;
        }
        return true;
    }
};

struct Scope : Term{
    Scope*    parrent;
    //CodeBlock body;
    //std::unordered_map<long,Type    *> funcs;
    //std::unordered_map<long,Function*> funcs;
    //std::unordered_map<long,Variable*> funcs;
    std::unordered_map<long,Term*> terms;

    Term* find_name( int i0, int n, long h )const{
        //long h = str_hash( n, source_code+i0, nameHash );
        auto got = terms.find(h);
        if(got != terms.end()){
            const Term& t = *got->second;
            if( t.equalName(i0, n) ) return (Term*)&t;
            printf( "ERROR Same hash %li, different name %.*s | %.*s \n", h, t.nsrc, source_code+t.isrc, i0, source_code+n );
        }
        return 0;
    }
    Term* find_name( int i0, int n )const{
        long h = str_hash( n, source_code+i0, nameHash );
        return find_name( i0, n, h );
    }

    Term* new_term( Term* t ){
        long h = str_hash( t->nsrc, source_code+t->isrc, nameHash );
        Term* found = find_name( t->isrc, t->nsrc, h );
        if(found){
            printf( "ERROR Term %.*s alrady defined \n", t->nsrc, source_code+t->isrc );
            return found;
        }
        terms.insert({h,t});
        return 0;
    }

};

struct Type : Term{
    Type* parrent=0;
    int   nbyte;
};

struct Variable : Term{
    void* data=0;
    int   nbyte;
    int   flags;

    //Variable(){}
};

struct Argument{
    Variable* var=0;
    int       flags;
};

struct Function : Term{
    std::vector<Argument> args;
    //CodeBlock body;
    Scope scope;
};

struct Operation{
    Function   func;
    std::vector<Variable*> operands;
    int        flags;
};

struct CodeBlock{
    std::vector<Operation> ops;
};

// ======== Main Program

class Program{

    int ich=0;
    long idCur        = 0; /// current ID
    //std::vector<Symbol>  symbols;
    Scope  scopeTop;
    Scope* scopeCur=&scopeTop;

    TokKind     tokKind = TokKind::none;
    Function*   curFunc=0;
    Type*       curType=0;
    Operation*  curOp=0;
    Argument*   curArg=0;

    int nextID(){
        return idCur++;
    }

    void parse( int nch ){
        int ichmax=ich+nch;
        for(ich=0;ich<ichmax;ich++){
            char ch = source_code[ich];
            if      (ch==ch_tokSep){ // token finished

            }else if(ch==ch_opSep ){ // operation finished
                switch(tokKind){
                    case TokKind::type:    break;
                    case TokKind::func:    break;
                    case TokKind::var:     break;
                    case TokKind::arg:     break;
                    case TokKind::literal: break;
                }
                tokKind = TokKind::none;
            }else if(ch==ch_def   ){ // new Term definition
                ich++;
                char ch = source_code[ich];
                if       (ch==ch_Func){ // new function
                }else if (ch==ch_Type){ // new variable
                }
            }
        }
    }

};


} // namespace ProgramGraph


#endif


