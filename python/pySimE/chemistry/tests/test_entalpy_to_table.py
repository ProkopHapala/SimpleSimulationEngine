#!/usr/bin/python

import re
import numpy as np
import sys
import ChemUtils as ch


'''
names = ( 'sucrose', 'potassium nitrate',   'potassium carbonate', 'carbon dioxide', 'water', 'nitrogen' )

print "============"
comps  = ch.names2comps( names, ch.substance_dict ) 
ch.printNiceMap( names, comps )

print "============"
masses = ch.comps2mass ( comps, ch.element_dict   )
ch.printNiceMap( names, masses )

print "============"
coefs  = ch.atomicBalance( names, comps, pivot=0 )
ch.printNiceMap( names, coefs )

print "============"
Hs = ch.takeSubstanceParam( names, 'enthalpy' )
ch.printNiceMap( names, Hs )

Htot = sum( coefs * np.array(Hs) )
print " enthalpy_change = ", Htot, " kJ/mol "

coefs_left  = -np.fmin( coefs, 0 );
coefs_right =  np.fmax( coefs, 0 )
Mtot        = sum( coefs_left * np.array( masses ) )
print " total mass = ", Mtot, " g/mol "

print " Energy density = ", Htot / Mtot, " MJ/kg "

#ch.atomicBalance( names, comps, pivot=0 )
'''



print
print 
print "============= sucrose + potassium nitrate "

print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'potassium nitrate',   'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'potassium nitrate',   'potassium hydroxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'potassium nitrate',   'potassium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )

print
print 
print "============= glycerol + potassium nitrate "
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'potassium nitrate',   'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'potassium nitrate',   'potassium hydroxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'potassium nitrate',   'potassium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )

print
print 
print "============= citric acid + potassium nitrate "
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'potassium nitrate',   'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'potassium nitrate',   'potassium hydroxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'potassium nitrate',   'potassium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )





print
print 
print "============= sucrose + calcium nitrate "
print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'calcium nitrate',   'calcium oxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'calcium nitrate',   'calcium oxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'sucrose', 'calcium nitrate',   'calcium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )

print
print 
print "============= glycerol + calcium nitrate "
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'calcium nitrate',   'calcium oxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'calcium nitrate',   'calcium oxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'glycerol', 'calcium nitrate',   'calcium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )

print
print 
print "============= citric acid + calcium nitrate "
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'calcium nitrate',   'calcium oxide', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'calcium nitrate',   'calcium oxide', 'carbon dioxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy(  [ 'citric acid', 'calcium nitrate',   'calcium carbonate', 'carbon dioxide', 'water', 'nitrogen' ] )



print
print 
print "============= aluminium + potassium permanganate "; ch.reactionEnthalpy(  [ 'aluminium', 'potassium permanganate', 'potassium aluminate', 'aluminium oxide', 'manganese' ] )
print "============= aluminium + potassium nitrate ";      ch.reactionEnthalpy(  [ 'aluminium', 'potassium nitrate', 'potassium aluminate', 'aluminium oxide', 'nitrogen' ] )
print "============= aluminium + ammonium  nitrate ";      ch.reactionEnthalpy(  [ 'aluminium', 'ammonium nitrate', 'aluminium oxide', 'water', 'nitrogen' ] )

print
print 
print "============= ammonium nitrate decompositon "; 
#print "============="; ch.reactionEnthalpy( [ 'ammonium nitrate', 'water', 'nitrous oxide' ] )
#print "============="; ch.reactionEnthalpy( [ 'ammonium nitrate', 'water', 'nitric oxide','nitrogen' ] )
#print "============="; ch.reactionEnthalpy( [ 'ammonium nitrate', 'water', 'nitrogen dioxide', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy( [ 'ammonium nitrate', 'water', 'nitrogen', 'oxygen' ] )
#print "============="; ch.reactionEnthalpy( [ 'ammonium nitrate', 'amonia', 'nitric acid' ] )



print
print 
print "============= ammonium chloride + ammonium nitrate"; 
print "============="; ch.reactionEnthalpy( [ 'amonium chloride', 'ammonium nitrate', 'water', 'nitrogen', 'hydrochloric acid' ] )

print
print 
print "============= carbon + ammonium nitrate "; 
print "============="; ch.reactionEnthalpy( [ 'carbon', 'ammonium nitrate', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy( [ 'carbon', 'ammonium nitrate', 'carbon dioxide',  'water', 'nitrogen' ] )

print
print 
print "============= sucrose + ammonium nitrate "; 
print "============="; ch.reactionEnthalpy( [ 'sucrose', 'ammonium nitrate', 'carbon monoxide', 'water', 'nitrogen' ] )
print "============="; ch.reactionEnthalpy( [ 'sucrose', 'ammonium nitrate', 'carbon dioxide',  'water', 'nitrogen' ] )




