#!/usr/bin/python

#!/usr/bin/python

'''
Variable
Expression =  Variable + Variable
Assignment(Expresssion)
'''

env = None
_float  = None
_int    = None
_string = None

#============================

class NameGen:

    def __init__(self,pre):
        #self.val = val
        self.pre = pre

    def __get__(self, instance, owner):
        return self.pre + instance.id

#============================

class Typ:

    def __init__(self,name):
        self.name = name
        self.id   = env.addType(self)

    def __str__(self):
        return self.name

#============================

class Operator():

    def __init__(self,name):
        self.name = name
        self.id = env.addOperator(self)

    def __str__(self):
        return self.name

#============================

class Constant:

    def __init__(self,val):
        self.val  = val
        if   isinstance(val, int  ): self.typ = env.types["int"]
        elif isinstance(val, float): self.typ = env.types["float"]
        elif isinstance(val, str  ): self.typ = env.types["string"]
        else: raise Exception( "ERROR Constant(x) not defined for type of x == ", type(val) )
        self.name = str(val)

class Var:

    def __init__(self,typ,name=None,val=None):
        self.val = val
        if   isinstance(typ, Typ):
            self.typ       = typ
        elif isinstance(typ, str):
            self.typ       = env.types[typ]
        self.name      = name
        #self.id        = env.get_next_var_id()
        self.id        = env.scope.add_var(self)
        if name is None :
            self.name = env.var_prefix + str(self.id)
            self.named = False
        elif isinstance(name, str):
            self.name  = name
            self.named = True
        else:
            raise Exception('ERROR in Var() : variable name must be string ' )

    def __add__(self,b): return Expression(env.operators['+'],self,b)
    def __sub__(self,b): return Expression(env.operators['-'],self,b)
    def __mul__(self,b): return Expression(env.operators['*'],self,b)
    def __div__(self,b): return Expression(env.operators['/'],self,b)
    def __pow__(self,b): return Expression(env.operators['**'],self,b)
    def __str__(self):   return self.name + " : "+self.typ.name

#============================

def checkArg(arg):
    if ( isinstance(arg, Var) or isinstance(arg, Expression) ): 
        return arg
    else: 
        return Constant(arg)

def checkOperator(op):
    if   isinstance(op, Operator): return op
    elif isinstance(op, Function): return op
    else:                          return Operator(op)

def checkStrUnfold(a):
    if isinstance(a, Expression):  return a.str_unfolded()
    else:                          return a.name

class Expression(Var):

    #name = NameGen("expr")

    def __init__(self,op, arg1, arg2=None, name=None ):
        self.op       = op
        if isinstance(self.op, Function):
            self.op       = op
            if isinstance(arg1, tuple):
                self.arg1     = arg1
            else:
                raise Exception( "ERROR Expression(op): Function argument-list must be tuple of types", type(arg1) )
            self.arg2     = None
        #elif op isinstance(self.op, Operator):
        else:
            self.op       = checkOperator( op )
            self.arg1     = checkArg( arg1 )
            self.arg2     = checkArg( arg2 )
        #else:
        #    raise Exception( "ERROR Expression(op): op is not Operator nor Function, but ", type(op) )

        #self.id       = env.get_next_var_id()
        self.id       = env.scope.add_expr(self)
        if name is None :
            self.name = env.expr_prefix + str(self.id)
            self.named = False
        else:
            self.name  = name
            self.named = True

#    def __add__(b):
#        if isinstance(b, Var) or isinstance(b, Expression):
#            return Expression(_add,self,b)
#        else:
#            raise Exception('ERROR in __add__( '+self+' , '+b+' ) : '+b+' is not variable nor expression ' )

    def str_assignment(self,folded=True):
        if folded:  return self.name+" : "+self.str_folded()
        else:       return self.name+" : "+self.str_folded()

    def str_folded(self):
        #print  "Expression.str_folded, types: ", type(self.arg1.name), type(self.operator), type(self.arg2.name)
        #print  "Expression.str_folded, types: ", type(self.arg1), type(self.operator), type(self.arg2)
        #print  "Expression.str_folded, types: ", self.arg1.__class__, self.operator.__class__, self.arg2.__class__
        #print  "    vals : ", self.arg1.name, str(self.operator), self.arg2.name
        if    isinstance(self.op, Operator):
            return "("+" ".join([self.arg1.name,str(self.op),self.arg2.name])+")"
        elif  isinstance(self.op, Function):
            return self.op.str_call(self.arg1)


    def str_unfolded(self):
        s1 = checkStrUnfold(self.arg1)
        s2 = checkStrUnfold(self.arg2)
        return "("+" ".join([ s1, str(self.operator), s2 ])+")"

    def __str__(self):
        if env.fold_expressions:  return self.str_folded  ()
        else:                     return self.str_unfolded()


#============================

class Scope(Expression):

    def __init__(self,parent,head=None):
        self.head = head
        self.parent = parent
        parent.childs.append(self)
        self.childs = []
        #self.variables   = parent.variables[:]
        #self.expressions = parent.expressions[:]  #   should we? expression are not named - we cannot refer to them
        #self.operators   = parent.operators.copy()
        #self.types       = parent.types.copy() 
        self.var_id0     = parent.var_id0  + len(parent.variables  )
        self.expr_id0    = parent.expr_id0 + len(parent.expressions)
        self.variables   = []
        self.expressions = []  #   should we? expression are not named - we cannot refer to them
        self.operators   = {}
        self.types       = {} 

    def add_var(self,a):
        self.variables.append(a)
        return self.var_id0  + len(self.variables)

    def add_expr(self,a):
        self.expressions.append(a)
        return self.expr_id0 + len(self.expressions)

    def str_head(self):
        s_indent = '    '*(self.level)
        if self.level==0:
            return s_indent+"GLOBAL_SCOPE {"
        else:
            return s_indent+str(self.head)

    def __str__(self):
        print self.level
        if env.fold_expressions:  return self.str_folded  ()
        else:                     return self.str_unfolded()

    def str_assignment(self,folded=True):
        if folded:  str_folded(self)
        else:       str_unfolded(self)

    def str_folded(self):
        s_indent = '    '*self.level
        return s_indent+str(head)

    def str_unfolded(self):
        s_indent = '    '*(self.level+1)
        ls = [ self.str_head() ]
        #s = [ "for %s in [%i..%i\%%i]{\n" %( self.i.name, self.i0.name, self.n.name, self.di.name) ]
        for e in self.expressions:
            if   isinstance(e, Scope):
                ls.append( e.str_unfolded() )
            elif isinstance(e, Expression):
                ls.append( s_indent + e.str_assignment() )
            else:
                ls.append( s_indent + str(e) )
        #print ls
        ls.append( '    '*(self.level) + "}")
        return "\n".join(ls)

#============================

class Range:

    def __init__(self,imax,i='i',i0=0,di=1):
        self.imax = checkArg(imax)
        self.i0   = checkArg(i0)
        self.di   = checkArg(di)
        if isinstance(i, str): self.i = Var(_int,i)
        else:                  self.i = checkArg(i0)
        self.bFirstIter = True

    def __iter__(self):
        self.bFirstIter = True
        return self

    def __str__(self):
        return "for %s in [%s..%s%%%s]{" %( self.i.name, self.i0.name, self.imax.name, self.di.name)

    def next(self):
        if self.bFirstIter:
            self.bFirstIter = False
            #print "RangeStart "
            env.beginScope(self)
            return self.i
        else:
            #print "RangeFinish "
            env.endScope()
            raise StopIteration

#============================

class Function:

    def __init__(self,name,argTypes=None,ret=None):
        self.name = name
        self.ret  = ret
        self.argTypes = argTypes
        self.args     = tuple( Var(t) for t in self.argTypes )

    def __iter__(self):
        self.bFirstIter = True
        return self

#    def __str__(self):
#        ls = [ self.name+"(" ]
#        for arg in self.args:
#            ls.append( arg.name )
#        return ",".join(ls) + ")"

    def __str__(self):
        ls = [ self.name+"(" ]
        #print args
        for i,(arg,typ) in enumerate(zip(self.args,self.argTypes)):
            if i>0 :
                ls.append(",")
            #if isinstance(arg, Var) or isinstance(i, Expression):
            #    ls.append( arg.name )
            #else:
            #    ls.append( str(arg) )
            #ls.append( str(arg)+":"+str(typ) )
            ls.append( str(arg) )
        ls.append( ")"+str(self.ret)+"{" )
        return "".join(ls)

    def next(self):
        if self.bFirstIter:
            self.bFirstIter = False
            #print "RangeStart "
            env.beginScope(self)
            return self.args
        else:
            #print "RangeFinish "
            env.endScope()
            raise StopIteration

    def __call__(self,*argv):
        return Expression(self, argv )

    def str_call(self,args):
        ls = [ self.name+"(" ]
        print args
        for i,arg in enumerate(args):
            if i>0 :
                ls.append(",")
            if isinstance(arg, Var) or isinstance(i, Expression):
                ls.append( arg.name )
            else:
                ls.append( str(arg) )
        return "".join(ls) + ")"

    def _return(self,e):
        pass

#============================

class Case:

    def __init__(self):
        pass

    def __eq__(self,val):
        env.scope.add_expr('==')
        print "Case .eq.",val

    def __ne__(self,val):
        print "Case .ne.",val

    def __gt__(self,val):
        print "Case .gt.",val

    def __ge__(self,val):
        print "Case .ge.",val

    def __lt__(self,val):
        print "Case .lt.",val

    def __le__(self,val):
        print "Case .le.",val

#  __lt__, __le__, __gt__, __ge__, __eq__ and __ne__

#============================

class Switch:

    def __init__(self,var):
        if isinstance(var, Var):
            self.var = var
        else:
            raise Exception( "ERROR Switch(var): var is not a Variable", type(var) )

    def __iter__(self):
        self.bFirstIter = True
        return self

    def next(self):
        if self.bFirstIter:
            self.bFirstIter = False
            #print "RangeStart "
            env.beginScope(self)
            return Case()
        else:
            #print "RangeFinish "
            env.endScope()
            raise StopIteration

    def __str__(self):
        return "Switch: "

#============================

class Environment(Scope):

    var_prefix  = "v"
    expr_prefix = "e" 

    def __init__(self):
        self.childs   = []   # from scope
        self.scopeSeq = []
        #self.next_var_id  = 0
        #self.next_expr_id = 0
        self.fold_expressions = False
        self.variables   = []
        self.expressions = []
        self.operators   = {}
        self.types       = {}
        self.rootScope   = self
        self.scope = self.rootScope
        self.scopeSeq.append(self.scope)
        self.level = 0
        self.var_id0     = 0
        self.expr_id0    = 0

    def initDefaults(self):
        global _int, _float, _string
        for s in ["int","string","float"]:
            Typ(s)
        _int    = self.types["int"]
        _float  = self.types["float"]
        _string = self.types["string"]
        for s in ["+","-","*","/","**",    "=="]:
            Operator(s)

    def addOperator(self, op ):
        self.scope.operators[op.name] = op
    
    def addType(self, typ ):
        self.scope.types[typ.name] = typ

    def beginScope(self,head):
        s          = Scope(self.scope,head)
        s.level    = 1 + self.scope.level
        self.scope = s
        #print self.scope.level, self.rootScope.level
        self.scopeSeq.append(s)
        self.expressions.append(s)

    def endScope(self):
        #print self.scope.level, self.rootScope.level
        if self.scope == self.rootScope:
            print "ERROR : cannot escape root scope "
        else:
            self.scope = self.scope.parent

#============================

if __name__ == "__main__":

    env = Environment()
    env.initDefaults()

    x = Var(_float,  "x")
    y = Var("float", "y")
    z = Var(_float,  "x")
    
    q = x ** 2
    
    a = x + y
    b = x - y
    
    f = a * b
    g = (f*2)*(a/b)
    
    rsum = Var(_float,'rsum',0)
    
    # loop over range
    for i in Range( 50 ):
        r2    = x*x

    # loop over range
    for i in Range( 20 ):
        rsum += x**2 + y**2 + z**2


    # define function
    f = Function("funkyFunction", (_float,_float,_float), _float )
    for x,y,z in f:
        h = (x/y + y/z)**3
        f._return( h )

    f(1,2,3)   # call the function

    for case in Switch(x):
        case == 5
        u = x - y*z
        case == 2
        u = z - x*z
        case == 3
        u = y - z*z


    print env.str_unfolded()

    exit()
    for i,v in enumerate(env.variables):
        print i,v
    
    for i,e in enumerate(env.expressions):
        print i,e
        #print i,e.str_assignment()













