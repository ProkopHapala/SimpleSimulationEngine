#!/usr/bin/python


thisProgram = None

class Type:

    def __init__(self,typ=None,val=None):
        self.typ = typ
        self.val = val

class Int(Type):


class Float(Type):


class Variable():

    def __init__(self,typ=None,val=None):
        self.typ = typ
        self.val = val

    def __set__(self, ins, v):
        if token.isinstance(Type):
            self.val
        elif token.isinstance(Variable):
            self.
        elif token.isinstance(Variable):

class Loop():

    def __init__(self,it):






class Program():

    def __init__(self):
        self.ops = []

    def __add__(self, op ):
        if op.isinstance(Variable):
            ops.append(op)




if __name__ == "__main__":
    thisProgram  = Program()
    i = Int()
    v = Vec(5,double,0)
    Loop(i,5, [
        q    = sqrt(i)  ,
        v[i] = q        ,
    ])
    import sys
    thisProgram.write(sys.stdout)






