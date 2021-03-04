
#ifndef  CommandTree_h
#define  CommandTree_h

#include "unordered_map"
#include "functional"

class CommandNode{ public:
    CommandNode* parrent=0;
    std::string info;
    std::function<void(void)> func;
    std::unordered_map<std::string,CommandNode*> leafs;

    CommandNode() = default;
    CommandNode( CommandNode* parrent_, std::string info_, std::function<void(void)> func_=[](){} ):parrent(parrent_),info(info_),func(func_){};

    CommandNode* addLeaf( std::string s, std::string info_, std::function<void(void)> func_=[](){} ){
        CommandNode* leaf = new CommandNode( this, info_, func_ );
        leafs[s] = leaf;
        return leaf;
    }

    char* leafsToStr(char* s){
        for( const auto& it : leafs ){
            s+=sprintf( s, "[%s] : %s \n",  it.first.c_str(), it.second->info.c_str() );
        }
        return s;
    }

};


class CommandTree{ public:
    CommandNode  root;
    CommandNode* cur =&root;

    void eval(std::string s){
        //printf( "DEBUG CommandTree::eval(%s)\n", s.c_str() );
        auto it = cur->leafs.find( s );
        if(it != cur->leafs.end()){
            cur = it->second;
            cur->func(); // call the lambda-function
        }
    }
    void goBack(){
        if(cur->parrent==0) return;
        cur=cur->parrent;
    }
    char* curInfo(char* s){
        return cur->leafsToStr(s);
    }

};

#endif
