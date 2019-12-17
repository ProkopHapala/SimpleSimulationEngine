
#ifndef CommandParser2_h
#define CommandParser2_h


/*

Simple scripting language desigmend for easy integration with strongly typed low-level language like C++
 - There are two disticnt types of functions
    - low-level callback ( which are type-checked )
    - hi-level macros which are type-less

## notation 1
  $out       = Func(in1,in2);
  $out1,out2 = Func(in1,in2);

## notation 2
  Func &out1 int1 &out2 int2 ;

## notation 3
  Func | out1 out2 | int1  int2 ;

*/


#include <vector>
#include <unordered_map>
#include <cstdio>
#include <cstring>

class VariableObject{ public:
    int   type;
    void* pointer;
    std::string name;
};

class FunctionObject{ public:
    int nArg     ;    //  =0;
    int retType  ;    //  =-1;
    int * argTypes;   // =0;
    std::string name;
};

inline bool isCapitalLetter (char c){ return (c>=65)&&(c<=90);  };
inline bool isSmallLetter   (char c){ return (c>=97)&&(c<=122); };
inline bool isDigit         (char c){ return (c>=48)&&(c<=57);  };
inline bool isTokenChar     (char c){ return isCapitalLetter(c)||isSmallLetter(c)||isDigit(c)||(c=='_');  };
inline bool isTokenSeparator(char c){ return (c==' ')||(c=='\t')||(c==',')||(c==';')||(c=='(')||(c==')'); };
inline bool isValueChar     (char c){ return isDigit(c)||(c=='.')||(c=='-')||(c=='+'); };

//==================

typedef std::vector<std::string> Strings;

class CommandParser2{ public:
    static constexpr int nArgMax = 16;
    static constexpr int nBufMax = 2048;
    static constexpr int maxTok  = 128;
    int       maxLines   = 10000;

    std::unordered_map<std::string,int> name2type;
    std::vector<std::string>            typeNames;

    std::unordered_map<std::string,int> name2funcion;
    std::vector<FunctionObject>         functions;

    std::unordered_map<std::string,int> name2variable;
    std::vector<VariableObject>         variables;

    bool bFuncNameChecking = true;

    bool bRet = false;
    std::string retName;
    std::string funcName;

    int   nArgs=0;
    void*       args      [nArgMax];
    double      argFloats [nArgMax];
    long        argInts   [nArgMax];
    std::string argStrings[nArgMax];
    int         argTypes  [nArgMax];

    char  buff    [nBufMax];
    int   iBuf = 0;

    char cComent  ='#';
    char cDefVar  ='$';
    char cDefFunc ='@';
    char cNewLine ='\n';
    //char cCmdEnd  ='\n';
    char cCmdEnd  =';';


    struct ArgType  { typedef enum { none=-1, Name=1, Float=2, Int=3,    String=4            } T; };
    struct TokenType{ typedef enum { none=-1, Name=1, Val=2,   String=3, Def=4,    Comment=5 } T; };
    struct NameType { typedef enum { none=-1, Func=1, Var=2                                  } T; };

    // ======= functions

    void printFunction(int i){
        printf( "%s %s(", typeNames[ functions[i].retType ].c_str(), functions[i].name.c_str() );
        for(int j=0; j<functions[i].nArg; j++){
            printf( " %s, ", typeNames[ functions[i].argTypes[j] ].c_str()  );
        }
        printf(")\n");
    };

    void printVariable(int i){
        printf( "variable[i].name=='%s' type '%s' &%i", i, variables[i].name, typeNames[variables[i].type], variables[i].pointer );
        switch(variables[i].type){
            case 0: printf( " %li \n", *(int*)  (variables[i].pointer) ); break;
            case 1: printf( " %lf \n", *(float*)(variables[i].pointer) ); break;
            case 2: printf( " %s  \n",  (char*) (variables[i].pointer) ); break;
            default: printf( "[Unknown] \n" ); break;
        }
    }

    void pushArgString(){
        buff[iBuf]='\0';
        argStrings[nArgs]=buff;
        args[nArgs]=&argStrings[nArgs];
        argTypes[nArgs]=2;
        //printf( "pushArg[%i] buff='%s' String '%s' '%s' \n", nArgs, buff, argStrings[nArgs].c_str(), ((std::string*)args[nArgs])->c_str() );
        nArgs++; iBuf=0;
    }

    int  pushArgFloat (){ buff[iBuf]='\0';
        int iret = sscanf( buff, "%lf", &argFloats[nArgs] );
        args[nArgs]=&argFloats[nArgs];
        argTypes[nArgs]=1;
        //printf( "pushArg[%i] buff='%s' Float %f' %f \n", nArgs, argFloats[nArgs], *(double*)args[nArgs] );
        if(iret!=1){printf("argument [%i] >>%s<< is not Float \n",nArgs,buff);};
        nArgs++; iBuf=0;
        return iret;
    }

    int  pushArgInt   (){ buff[iBuf]='\0';
        int iret = sscanf( buff, "%li", &argInts[nArgs] );
        args[nArgs]=&argInts   [nArgs];
        argTypes[nArgs]=0;
        //printf( "pushArg[%i] buff='%s' Int %li %li \n", nArgs, buff, argInts[nArgs], *(long*)args[nArgs] );
        if(iret!=1){printf("argument [%i] >>%s<< is not Int   \n",nArgs,buff);};
        nArgs++; iBuf=0;  return iret;
    }

    int pushArgName  (){ buff[iBuf]='\0';
        std::string vname = buff;
        auto got = name2variable.find(vname);
        if ( got == name2variable.end() ){ printf( "Varaible '%s' not known \n", vname.c_str() ); return -1; }
        else{
            VariableObject&  var = variables[got->second];
            //printf( "pushArgName '%s' variables[%i].name='%s' of type '%s' &%i  \n", vname.c_str(), got->second, var.name.c_str(), typeNames[var.type].c_str(), var.pointer );
            args[nArgs]     = var.pointer;
            argTypes[nArgs] = var.type;
        };
        nArgs++; iBuf=0;
        return 0;
    }

    int execFunction(){
        // check function exist
        //printf( "execFunction : '%s'[%i] -> '%s' \n", funcName.c_str(), nArgs, retName.c_str() );
        auto gotf = name2funcion.find(funcName);
        if  ( gotf == name2funcion.end() ){ printf( "Functions >>%s<< not known \n", funcName.c_str() ); return -1; }
        int ifunc = gotf->second;
        //printFunction(ifunc);
        FunctionObject& func = functions[ifunc];
        // check argument number and types
        if(func.nArg!=nArgs){ printf( "Functions >>%s<< should have %i Arguments (%i given) \n", funcName.c_str(), func.nArg, nArgs ); return -2; }
        for(int i=0; i<func.nArg; i++){
            if( func.argTypes[i] != argTypes[i] ){
                printf("Function >>%s<< argument[%i] type is >>%s<< (should be >>%s<<) \n", funcName.c_str(), i, typeNames[argTypes[i]].c_str(), typeNames[func.argTypes[i]].c_str() );
                return -(i+2);
            }
        }
        // check return value exist
        if(bRet){ if ( name2variable.find(retName) != name2variable.end() ){ printf( "Variable >>%s<< already defined \n", retName.c_str() ); return -3; } }
        // dispatch call to C++
        void* ret = callTable(ifunc,nArgs, args );
        // store return value
        if( bRet ){
            name2variable.insert({retName,variables.size()});
            variables.push_back({func.retType,ret,retName});
        }
        return 0;
    }

    inline void clearExecTemps(){ nArgs=0;bRet=false; }

    void execString( int n, char* line ){

        TokenType::T inWhat   = TokenType::none;
        NameType ::T nameType = NameType ::none;
        bool isFunc    =0; // is this name of function ?
        bool isInt     =0; // is this value integer ?

        for(int i=0; i<n; i++ ){
            char c = line[i];
            //printf( "%i '%c' %i ", i, c, inWhat );
            if( c=='\0' ) return;
            switch(inWhat){
                case TokenType::Name:
                    //printf("in-Name " );
                    if(isTokenChar(c)){ buff[iBuf++]=c; }
                    else{
                        switch( nameType ){
                            case NameType::Func: buff[iBuf]='\0'; funcName=buff;  iBuf=0; break;
                            case NameType::Var:  pushArgName(); break;

                        }
                        inWhat = TokenType::none;
                        //if( c=='"' ){ inWhat = TokenType::string; }else{ inWhat = TokenType::string; }
                    };
                    break;
                case TokenType::Def:
                    //printf("in-Def " );
                    if(isTokenChar(c)){ buff[iBuf++]=c; }
                    else{
                        switch( nameType ){
                            //case NameType::func: break;
                            case NameType::Var: buff[iBuf]='\0'; retName=buff; bRet=true; iBuf=0;  break;
                        }
                        inWhat = TokenType::none;
                        //if( c=='"' ){ inWhat = TokenType::string; }else{ inWhat = TokenType::string; }
                    };
                    break;
                case TokenType::Val:
                    //printf("in-Val " );
                    if( isValueChar(c) ){ buff[iBuf++]=c; if(c=='.')isInt=false; }
                    else{
                        if(isInt){ pushArgInt(); }else{ pushArgFloat(); };
                        inWhat = TokenType::none;
                    };
                    break;
                case TokenType::String:
                    //printf("in-String " );
                    if(c=='"'){ pushArgString(); inWhat=TokenType::none; }else{ buff[iBuf++]=c; };
                    break;
                case TokenType::Comment:
                    //printf("in-Comment " );
                    if( (c==cComent)||(c==cNewLine) ){ inWhat=TokenType::none; };
                    break;
                default:
                    if     ( c==cCmdEnd         ){ execFunction(); clearExecTemps(); }
                    else if( c==cComent         ){ inWhat=TokenType::Comment; } // break For cycle
                    else if( c=='"'             ){ inWhat=TokenType::String;  }
                    else if( isValueChar(c)     ){ inWhat=TokenType::Val;  isInt=true;              buff[iBuf++]=c; }
                    else if( isCapitalLetter(c) ){ inWhat=TokenType::Name; nameType=NameType::Func; buff[iBuf++]=c; }
                    else if( isSmallLetter(c)   ){ inWhat=TokenType::Name; nameType=NameType::Var;  buff[iBuf++]=c; }
                    else if( c==cDefVar         ){ inWhat=TokenType::Def;  nameType=NameType::Var;  } // variable definition
                    //if( c==cDefFunc ){ inWhat = TokenType::name; nameType=NameType::var; } // variable definition
            }
            //printf( "\n" );
        }
    }

    int registerType     (std::string str){
        auto got = name2type.find(str);
        if ( got == name2type.end() ){
            int i=typeNames.size();
            name2type.insert( {str,i} );
            typeNames.push_back(str);
            return 1;
        }else{
            printf( "Type >>%s<< already known => IGNORED \n", str.c_str() );
            return -1;
        };
    }

    // https://stackoverflow.com/questions/9626722/c-string-array-initialization
    int registerFunction (std::string fname, std::string retType, Strings const& argTypes ){
        int iret;
        auto gotf = name2funcion.find(fname);
        if ( gotf != name2funcion.end() ){ printf( "Function >>%s<< alread known (%i-th) \n", fname.c_str(), gotf->second ); return -1; }
        auto got = name2type.find(retType);
        if ( got == name2type.end() ){ printf( "Function >>%s<< return type >>%s<< not known \n", fname.c_str(), retType.c_str() );     return -2; }
        else                         { iret = got->second; };
        int* argTypeIs = new int[argTypes.size()];
        int i=0;
        for( std::string argType: argTypes ){
            got = name2type.find(argType);
            if ( got == name2type.end() ){ printf( "Function >>%s<< argument [%i] type >>%s<< not known \n", fname.c_str(), i, argType.c_str() ); return -(i+2); }
            else                         { argTypeIs[i] = got->second; };
            i++;
        }
        name2funcion.insert( {fname,functions.size()} );
        functions.push_back( {argTypes.size(),iret,argTypeIs,fname} );
        printFunction(functions.size()-1);
        return 0;
    }

    inline int checkFunc(int i, const char * fname){
        //printf( "checkFunc functions[%i].name=='%s' | fname='%s' \n", i, functions[i].name.c_str(), fname );
        if(bFuncNameChecking){
            int icmp = strcmp( fname, functions[i].name.c_str() );
            if(icmp!=0) printf( "ERROR functions[%i].name=='%s' not '%s' \n", i, functions[i].name.c_str(), fname );
            return icmp;
        }
        return -1;
    }

    virtual void* callTable( int ifunc, int nArgs, void** args )=0;

    CommandParser2(){
        name2type.insert( {"void",-1} );
        //name2type.insert( {"void",-1} );
        //name2type.insert( {"void",-1} );
        registerType( "Int"    );
        registerType( "Float"  );
        registerType( "String" );
    }

};

#endif


