
#include "StructParser.h"
#include "AST.h"
#include "AST2Lang.h"

#include "binarySwitch.h"
#include "testUtils.h"


#include <cstring>

using namespace Lang;
using namespace parsing;

StructParser parser;
AST ast;


#include "TreeParser.h"
TreeParser tparser;


//char * test_str = "Name1{ x1; [10,15,5.5,-18.9]; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

char * test_str = "Name1{ x1; x2; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";

//template<class T,int N>
//struct VecNT{ T arr[N]; };


//template<class T,int N>
//T sumNT(T* xs){ T sum=0; for(int i=0;i<N;i++){sum+=xs[i];}; return sum; };


int main(){




    tparser.checkTables();

    //char * test_str = "Name1{ x1; x2; x3; Name2{ x3; x4; };  Name3{ x5; x6; };  };";
    //char * test_str = " Name1 + Name2 + Mame3 ";
    //char * test_str = " AAA + BBB; + CCCasd += DDD158 ; XXasdX , Y16YY";

    //char * test_str = " AA , BB ; CC , DD ";

    //char * test_str = " AA*H , U+BB ; CC+B , DD, XX ; MM ";
    //char * test_str = " AA*H , U+BB ";

    char * test_str = " AA + ( BB + CC ) ";

    //char * test_str = " AAA + BBB + CCC ";
    //char * test_str = " AAA BBB CCC ";

    tparser.parseString( strlen(test_str), test_str );
    printf("\n");
    printf("%s\n", tparser.str);
    tparser.evalAllLevels();


    //printf("\033[0;31m");

    tparser.printItems();
    tparser.printLevelColor();

    exit(0);





    //printf( " sizeof(std::vector<int> %i \n", sizeof(std::vector<int>) );
    //printf( " sizeof(std::unordered_map<int,int>) %i \n", sizeof(std::unordered_map<int,int>) );
    //printf( " sizeof(std::unordered_map<std::string,int>) %i \n", sizeof(std::unordered_map<std::string,int>) );

    //using Vec5d = VecNT<double,5>;

    //using sum5d = sumNT<double,5>;
    //using sum5d(double* xs) = sumNT<double,5>(T* xs);

    //using binarySwitch_(const T c,const T* levels) = binarySwitchT<int,3,9>(const T c,const T* levels);
    //using binarySwitch_ = binarySwitchT<int,3,9>;

    //  " !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    char ascii[] = "                                                                                                                                              ";
    for(int i=32;i<128;i++){
        ascii[i]=i;
    }
    printf("\"%s\"\n", ascii );

    printf("0 in ascii %i %c \n", (int)'0', '0' );

    //int levels[]{2,5,8,16,20,32,45,55, 69,80,0x7FF,0x7FF,0x7FF,0x7FF,0x7FF,0x7FF};

    const int mlev = 3;
    const int nlev = 2<<mlev;
    int levels[nlev];
    int lev = 3;
    for(int i=0;i<nlev; i++){ levels[i] = lev; lev+=rand()%5; };
    //int c = rand()%;
    /*
    int c = 20;
    printf( "c : %i \n", c );
    //size_t i = binarySwitch<int,2,8>( c, levels );4
    //int i = binarySwitch( c, 3, 9, levels);
    int i = binarySwitchT<int,4,9>( c, levels);
    if(i>0){ printf( "[%i]:  %i == %i      \n", i, c, levels[i] ); }
    else   { printf( "[%i]:  %i not found  \n", i, c            ); }
    */



    //const char s    [] = "for(int i=0;i<5;i++;){a[i]=sqrt(b[i]);};";
    /*
    const char cases[] = "()[]{}\255\255";
    const char* pc = s;
    while((*pc)!=0){
        int icase = binarySwitchT<char,2,6>( *pc, cases );
        if( icase>0 )printf( "s[%i] case[%i] found %c \n", pc-s, icase, *pc );
        pc++;
    }
    */

    // "123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

    //const char   from[] = "({[ )}] ,; .+-*/%^~&|?<>=:!#$@0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_\t\n";
    //const char   to[] = "((( ))) ,, +++++++++++++++++++1111111111AAAAAAAAAAAAAAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaaaaaaaaaaa  ";
    //const char   to[] = "((( ))) ,, +++++++++++++++++++aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa  ";
    //const char     to[] = "((( ))) ,, +++++++++++++++++++aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa  ";

    //char table[255];
    //bakeSwitchTable<char,char>( 255, strlen(from), from, to, table, 0, 0, ' ' );


    /*
    int nchs = strlen(s);
    char s_[nchs];
    for(int i=0;i<nchs;i++){ s_[i]=table[s[i]]; }
    */


    /*
    printf("============\n" );
    char s_[255];
    for(int i=0;i<255;i++){ s_[i]=' ';       };
    for(int i=33;i<127;i++){ s_[i]=i;        }; printf("%s\n", s_ );
    for(int i=33;i<127;i++){ s_[i]=table[i]; }; printf("%s\n", s_ );
    printf("============\n" );
    exit(0);


    int n = 1000;
    int m = 1000;
    int* arr = new int[n];
    for(int i=0; i<n; i++){ arr[i]=rand()%0x7F; };
    //SPEED_TEST_FUNC_NM( "template    ", sqrt(16), n, m );
    //SPEED_TEST_FUNC_NM( "template    ", (binarySwitchT<int,3,9>( c,      levels)) , n, m );
    //SPEED_TEST_FUNC_NM( "no-template ",  binarySwitch ( c, 3,9, levels)            , n, m );

    SPEED_TEST_FUNC_NM( "template    ", (binarySwitchT<int,mlev,nlev>( arr[i],            levels) ) , n, m );
    SPEED_TEST_FUNC_NM( "no-template ", (binarySwitch                ( arr[i], mlev,nlev, levels) ) , n, m );

    //printf( "[%i]:  %i < %i < %i \n", i, levels[i], c, levels[i+1] );

    //printf( "%i \n", i );
    exit(0);
    */

    /*
    printf( " %s\n", test_str );
    printf( " ==== \n" );
    parser.parseString( strlen(test_str), test_str );
    //parser.parse( 0, 0, -1 );
    parser.printItemStruct();
    */

}
