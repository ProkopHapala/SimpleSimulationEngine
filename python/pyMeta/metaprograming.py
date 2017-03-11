#!/usr/bin/python

import os
import datetime
import re

lib_pre = 'lib'

re_curly = re.compile( '(?:|{})', re.DOTALL)
#re_currly = re.compile( '(?:|{})', re.DOTALL)


# see   https://docs.scipy.org/doc/numpy/user/basics.types.html
array_types = {
 "float"    : ( "float32"  , "array1f"   ),
 "double"   : ( "float64"  , "array1d"   ),
 "int"      : ( "int32"    , "array1i"   ),
 "bool"     : ( "int"      , "array1b"   ), 
 "unit8_t"  : ( "uint8"    , "array1u8"  ),
 "unit16_t" : ( "uint16"   , "array1u16" ),
 "unit32_t" : ( "uint32"   , "array1u32" ),
 "unit64_t" : ( "uint64"   , "array1u64" ),
 "char"     : (  None      , "c_char_p"  ),
}

# see : https://docs.python.org/2/library/ctypes.html
prim_types = {
"float"    : "c_float", 
"double"   : "c_double", 

"bool"     : "c_bool", 		 	
"char"     : "c_char",	 	

"int" 	   : "c_int",
"short"    : "c_short",
"long"     : "c_ulong",
 
"uint8_t"  : "c_ubyte",
"uint16_t" : "c_ushort",
"uint32_t" : "c_uint",
"uint64_t" : "c_ulong",
	 	
#"wchar_t"            : "c_wchar"
#"unsigned char"      : "c_ubyte" 	
#"short"              : "c_short" 
#"unsigned short"     : "c_ushort" 
#"unsigned int" 	  : "c_uint"
#"long"               : "c_long"
#"unsigned long"      : "c_ulong" 
#"__int64"            : "c_longlong"
#"long long"          : "c_longlong" 
#"unsigned __int64"   : "c_ulonglong"
#"unsigned long long" : "c_ulonglong"
#"long double" 	      : "c_longdouble" 
}

old_array_types = set()
new_array_types = set()

'''
# http://stackoverflow.com/questions/29974929/match-open-and-close-brackets-in-garbled-string

def tokenize(source):
    start = 0
    end = len(source)
    while start < end:
        match = pattern.match(source, start)
        if match:
            yield match.group(0)
        else:
            raise ValueError('Invalid syntax at character %d' % start)
        start = match.end()

def findCloseBrace( lines, iline0, ichr0, q='{', p='}' ):
    brlevel = 1
    iline   = iline0 
    ichr    = ichr0
    while brlevel>0:
        print lines[iline],
        tokens = re_curly.match( lines[iline] ).group(0)
        print tokens
        iline+=1
        if(iline>iline0+100):
            break
'''


def parse_cpp_header( iline0, ichr, lines ):
    args = []
    nextl=True
    iline = iline0
    iend  = iline0+100
    ichr  = lines[iline].find('(',ichr)+1
    while(iline<iend):
        line = lines[iline][ichr:];  ichr=0 
        ibr  = line.find(')')
        if(ibr>-1):
            line   = line[:ibr]
            iend   = 0
        else:
            iline += 1
        print "line %i: >>%s<<" %(iline,line)
        args += line.strip().split(',')
    #args = filter("", args)
    args = [x for x in args if x!=""]
    print "args : ", args
    return args,iline,ibr
    
def parse_cpp_func( iline0, ichr, fun_name, lines ):
    ret=lines[iline0][:ichr].split()[-1]
    args,iline,ichr = parse_cpp_header( iline0, ichr+len(fun_name), lines )
    func_obj = [fun_name,args,ret,lines[iline0:iline+1]]
    return iline, func_obj

def ctype2py( cstr ):
    global new_array_types
    ichr  = cstr.find('*')
    if(ichr>-1):
        for key,val in array_types.iteritems():
            if key in cstr:
                ichr+=1
                if not (val in old_array_types):
                    new_array_types.add(val)
                    #print "adding arraytype : ", val,new_array_types 
                return (cstr[ichr:].strip(),cstr[:ichr].strip(),val[1])
    else:
        for key,val in prim_types.iteritems():
            ichr = cstr.find(key)
            if ichr>-1:
                ichr+=len(key)
                return (cstr[ichr:].strip(),cstr[:ichr].strip(),val)
    print "ERROR ctype2py:   ", cstr
    #return "TYPE_NOT_FOUND"
    return None
    
def write_python_interface(func_obj, wraper=True,  cppcomment=True ):
    funname= func_obj[0]
    #print "args : ", func_obj[1]
    args     = [ ctype2py( arg ) for arg in func_obj[1] ]
    ret      = prim_types.get(func_obj[2])
    print  funname,args
    s = "" 
    if cppcomment:
        s += "".join( [ "#"+s for s in func_obj[3] ] )                                                        
    s +=  "%s.%s.argtypes = [%s]\n" %(lib_pre,funname,",".join([arg[2] for arg in args]) )
    s +=  "%s.%s.restype  = [%s]\n" %(lib_pre,funname,ret)                                
    if( wraper ):
        argnames = ",".join([arg[0] for arg in args]) 
        if ret is not None:
            retstr="return"
        else:
            retstr=""
        s +=(("def %s( %s ):\n" %(funname,argnames))+
            ("   %s %s.%s( %s )" %(retstr,lib_pre,funname,argnames))
        )
    return s

def write_arraytypes( ):
    return "".join([ "%s = np.ctypeslib.ndpointer(dtype=np.%s, ndim=1, flags='CONTIGUOUS')\n" %(typ[1],typ[0]) for typ in new_array_types ])
            
def check_existing_tokens( func_names, py_name ):
    global old_array_types
    func_names_ = [ "lib.%s.argtypes" %s for s in func_names ]
    with open(py_name,'r') as py_file:
        for line in py_file:
            for fun_name in func_names_:
                if fun_name in line:
                    func_names.remove(fun_name)
            if( "np.ctypeslib.ndpointer" in line ):
                for typ in array_types.values():
                    if (typ[1] in line):
                        old_array_types.add(typ)
                        break

def generate_interface( cpp_name, py_name, func_names ):
    if( os.path.isfile(py_name) ):
        check_existing_tokens( func_names, py_name )
        print "old_array_types", old_array_types
        print "func_names", func_names
        py_file  = open(py_name, 'a')
        py_file.write("#========= auto update : %s\n" %str(datetime.datetime.now()) )
    else:
        py_file  = open(py_name, 'w')
    with open(cpp_name) as f:
        cpp_lines = f.readlines()
    iline = 0
    nc = len(cpp_lines)
    new_funcs_interfaces = []
    while (iline < nc):
        cline =  cpp_lines[iline]
        #print iline, cline,
        for fun_name in func_names:
            ichr = cline.find(fun_name)
            if ichr > -1:
                #print ">>>>>>> FOUND : ",fun_name 
                iline,func_obj = parse_cpp_func( iline, ichr, fun_name, cpp_lines ) 
                new_funcs_interfaces.append( write_python_interface(func_obj) )
                func_names.remove(fun_name)
                break
        iline+=1 
    print "old_array_types", old_array_types
    print "new_array_types", new_array_types
    py_file.write( write_arraytypes( ) + "\n" ) 
    for func in new_funcs_interfaces:
        py_file.write( func + "\n" )
        py_file.write("\n")  
    py_file.close

    
