
#include "StructParser.h"
#include "AST.h"
#include "AST2Lang.h"

#include <cstring>

using namespace Lang;

StructParser parser;
AST ast;


//char * test_str = "Name1{ x1; [10,15,5.5,-18.9]; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

char * test_str = "Name1{ x1; x2; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

int main(){
    printf( " sizeof(std::vector<int> %i \n", sizeof(std::vector<int>) );
    printf( " sizeof(std::unordered_map<int,int>) %i \n", sizeof(std::unordered_map<int,int>) );
    printf( " sizeof(std::unordered_map<std::string,int>) %i \n", sizeof(std::unordered_map<std::string,int>) );

    printf( " %s\n", test_str );
    printf( " ==== \n" );
    parser.parseString( strlen(test_str), test_str );
    //parser.parse( 0, 0, -1 );
    parser.printItemStruct();
}
