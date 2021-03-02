
#ifndef  repl_h
#define  repl_h

#include <any>
#include <functional>
#include <string>
#include <unordered_map>

namespace REPL{ // https://en.wikipedia.org/wiki/Read_eval_print_loop

void init(){
	// Non-Blocking terminal input : see https://stackoverflow.com/questions/6055702/using-fgets-as-non-blocking-function-c
    int fd     = fileno(stdin);
    int flags  = fcntl(fd, F_GETFL, 0);
    flags |= O_NONBLOCK;
    fcntl(fd, F_SETFL, flags);
}


class Interpreter{ public:
    // https://stackoverflow.com/questions/61969316/is-it-possible-to-put-lambda-expressions-into-a-map-or-list-in-c
    // https://en.cppreference.com/w/cpp/utility/any
    std::unordered_map< std::string, std::function<void()>> functions;
    //std::unordered_map< std::string, std::function<std::any(std::any)>> functions;

    void eval(){
        char inputLine[1024];
        char* s = fgets(inputLine,1024,stdin);

        if(s!=0){
            auto it = functions.find(s);
            if ( it != functions.end() ){
                (it->second)();
            }else{
                printf( "REPL: %s No such command \n", s );
                printf( "REPL: knows: \n", s );
                for( auto it : functions ){ printf( "functions[] %s \n", it.first.c_str() );  };
            }
        }
    }

};

}; //  namespace REPL

#endif
