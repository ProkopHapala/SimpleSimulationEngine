
#include "StructParser.h"

#include "binarySwitch.h"
#include "testUtils.h"

#include <cstring>

using namespace parsing;

StructParser parser;


char * test_str = "Name1{ x1; x2; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

int main(){


    printf( " %s\n", test_str );
    printf( " ==== \n" );
    parser.parseString( strlen(test_str), test_str );
    //parser.parse( 0, 0, -1 );
    parser.printItemStruct();
}
