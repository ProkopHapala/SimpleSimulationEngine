#!/usr/bin/python

import os 

recompile = True 
ext       ='.so'

def compile_lib( name,
		FFLAGS = "-std=c++11 -Og -g -Wall",
		LFLAGS = "-I/usr/local/include/SDL2 -lSDL2",
		path   = None,
		clean  = True,
	):
	print " COMPILATION OF : "+name
	if path is not None:
		dir_bak = os.getcwd()
		os.chdir( path);
	print os.getcwd()
	if clean:
		try:
			os.remove( 'lib'+name+ext )
			#os.remove( name+".o" ) 
		except:
			pass 
	os.system("g++ "+FFLAGS+" -c -fPIC "+name+".cpp -o "+name+".o "+LFLAGS )
	os.system("g++ "+FFLAGS+" -shared -Wl,-soname,lib"+name+ext+" -o lib"+name+ext+" "+name+".o "+LFLAGS)
	if path is not None:
		os.chdir( dir_bak )
