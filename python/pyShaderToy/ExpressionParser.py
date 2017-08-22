#!/usr/bin/python

import sys
sys.stdout.write('.')

OPERATORS = { 
'+' : 0, 
'-' : 1, 
'*' : 2, 
'/' : 3, 
'^' : 4, 
}

def tokenTree( s, i0, imax, level, cOPEN="(",cCLOSE=")" ):
    tabs = " "*4*level
    toks = []
    expr = [] 
    i    = i0
    oi   = i0 
    nonEmpty = False;
    curr_op  = None

    s_debug = ""

    def finishSimpleTok(c):
        global nonEmpty
        global s_debug
        tok = (curr_op,(oi,i))  #    ( operator, (istart,iend) )
        toks.append( tok ); 
        nonEmpty = False;
        #print tabs+"(%i,%i)  :  " %(oi,i), s[oi:i]
        #s_debug = tabs+"(%i,%i)  :  " %(oi,i), s[oi:i]
        
    while i<imax:
        c = s[i]
        #s_debug = c+"|    "
        if c == cOPEN:
            iend,t  = tokenTree(s, i+1, imax, level+1, cOPEN, cCLOSE )
            #print tabs+"(%i,%i)  :  " %(oi,iend), s[oi:iend]
            #s_debug =  tabs+"(%i,%i)  :  " %(oi,iend), s[oi:iend]
            tok = (curr_op, (oi,i,iend), t )
            toks.append( tok );      #    ( operator, (istart,iend), function,arguments )
            i  = iend
            oi = i
            curr_op  = None
            nonEmpty = False
        elif c == cCLOSE:
            if nonEmpty: finishSimpleTok(c)
            return i+1, toks
        elif c in OPERATORS:
            if nonEmpty: finishSimpleTok(c)
            curr_op = c
            oi      = i+1
        elif (c ==" ") or (c =="\t") or (c =="\n"):
            pass
        else:
            nonEmpty=True
        #print s_debug
        i+=1
    return imax, toks

def walkTree_1( s, tree, level, fItem, fPre, fPost ):
    for t in tree:
        tlen = len(t)
        if tlen == 2:
            fItem(s,t,level)
        if tlen == 3:
            fPre(s,t,level)
            walkTree_1( s, t[2], level+1, fItem, fPre, fPost )
            fPost(s,t,level)

def walkTree_2( s, tree, level, fItem, fPre, fPost ):
    for t in tree:
        tlen = len(t)
        op = t[0]
        it = t[1]
        st = s[it[0]:it[1]]
        if tlen == 2:
            fItem(op,st,level)
        if tlen == 3:
            fPre(op,st,level)
            walkTree_2( s, t[2], level+1, fItem, fPre, fPost )
            fPost(op,st,level)

def tokenTree2string( s, tree, level, cOPEN="(",cCLOSE=")", indent=False ):
    so = ""
    if indent: 
        tabs = " "*4*level
        cCLOSE_ = tabs+("%s\n" %cCLOSE )
    else: cCLOSE_=cCLOSE;  tabs="";
    
    def putTok(it,op):
        st = s[it[0]:it[1]].strip()
        if op is None:
            st=" "+st
        else:
            st=op+st
        return st
        
    for t in tree:
        #print t
        tlen = len(t)
        if tlen == 2:
            st = putTok(t[1],t[0])
            if indent: st = tabs+("%s\n" %st )
            so += st
        if tlen == 3:
            st = putTok(t[1],t[0])+cOPEN
            if indent: st =  "%s%s\n" %(tabs,st)
            so += st + tokenTree2string( s, t[2], level+1, cOPEN,cCLOSE, indent )+cCLOSE_
    return so

if __name__ == "__main__":
    s = "sin( x + sin( (y * 6) - z ) ) + x*cos( (u * w) - c )"
    
    print s
    print "================"
    iend,tree = tokenTree( s, 0, len(s), 0 )
    print "================"
    print tree
    print "================"
    ss =  tokenTree2string( s, tree, 0, indent=True )
    print "================"
    print ss
    
    
    
    
    
    
