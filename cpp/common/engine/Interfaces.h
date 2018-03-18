
#ifndef Interfaces_h
#define Interfaces_h

class AnyControler{ public:
    //virtual void control(void* obj, double dt)=0;
    virtual void update(double dt)=0;
};

#endif

