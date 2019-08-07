
#ifndef  str2tree_h
#define  str2tree_h

class Str2Tree{ public:
    static const int maxLevel = 8;
    int level =0;
    int cur [8];
    int goal[8];
    char copen='{',cclose='}',csep=';';
    int ich    =0,nch=0;
    int nchMax = 1000000;
    const char* str=0;

    void setStr(const char* str_){str=str_;ich=0;nch=0;level=0;}

    //inline int branch(char c){}

    int step( char* stmp ){
        ich+=nch;
        nch=0;
        int ret;
        while(true){
            char c = str[ich+nch];
            nch++;
            //printf( "%c \n", c );
            if      ( c==csep   ){ cur[level]++;          ret= 0; break; }
            else if ( c==copen  ){ level++; cur[level]=0; ret=-1; break; }
            else if ( c==cclose ){ level--;               ret= 1; break; }
            //else if (stmp){stmp}
            else if(stmp){ *stmp=c; stmp++; }
            //if(ich>nchMax) break;
        }
        if(stmp)*stmp='\0';
        return ret;
    }

    void seekN(int n){ for(int i=0; i<n; i++)step(0); }

    int checkGoal(int nlevel){
        for(int i=0; i<nlevel;i++){
            int ci=cur[i];
            int gi=goal[i];
            if     (ci>gi ){ return  i; }
            else if(ci!=gi){ return -i; }
            //if     (cur[i]>goal[i] ){ return -i; }
        }
        return 0;
    }

    int seekGoal(int nlevel){
        int ret;
        while(true){
            step(0);
            ret = checkGoal(nlevel);
            if(ret>=0) return ret;
        }
        return ret;
    }

    /*
    void seek(int level_, int cur_){
        do{ step(0); }while( (level!=level_)||(cur[level]!=cur_) );
        //do{ step(0); if((level==level_)&&(cur[level]==cur_))break; };
        //while( (level!=level_)||(cur[level]!=cur_) ){  };
    }
    */
};

#endif
