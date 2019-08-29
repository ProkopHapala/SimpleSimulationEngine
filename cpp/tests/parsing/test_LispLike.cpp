#include "StructParser.h"
//#include "LispParser.h"

#include <cstring>

using namespace parsing;


StructParser parser;

//LispParser parser;



//char * test_str = "Name1{ x1; [10,15,5.5,-18.9]; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

//char * test_str = "Name1{ x1; x2; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";
char * test_str = "N1{ x1; x2; x3; N2{ x3; x4; N3{ x7; x8; }; };  N4{ x5; x6; };  };";

int main(){
    printf( " %s\n", test_str );
    printf( " ==== \n" );
    parser.parseString( strlen(test_str), test_str );
    //parser.parse( 0, 0, -1 );
    parser.printItemStruct();
}
