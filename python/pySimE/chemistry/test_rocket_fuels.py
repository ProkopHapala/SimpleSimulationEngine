#!/usr/bin/python

import re
import numpy as np
import sys
import ChemUtils as ch





#        Fuel           Oxidizer               Producs      

reactions=[
( 'sucrose',          'potassium nitrate',  [ 'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] ),
( 'sucrose',          'potassium nitrate',  [ 'potassium hydroxide', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'sucrose',          'potassium nitrate',  [ 'potassium carbonate', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'glycerol',         'potassium nitrate',  [ 'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] ),
( 'glycerol',         'potassium nitrate',  [ 'potassium hydroxide', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'glycerol',         'potassium nitrate',  [ 'potassium carbonate', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'citric acid',      'potassium nitrate',  [ 'potassium hydroxide', 'carbon monoxide', 'water', 'nitrogen' ] ),
( 'citric acid',      'potassium nitrate',  [ 'potassium hydroxide', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'citric acid',      'potassium nitrate',  [ 'potassium carbonate', 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'sucrose',          'calcium nitrate',    [ 'calcium oxide',       'carbon monoxide', 'water', 'nitrogen' ] ),
( 'sucrose',          'calcium nitrate',    [ 'calcium oxide',       'carbon dioxide',  'water', 'nitrogen' ] ),
( 'sucrose',          'calcium nitrate',    [ 'calcium carbonate',   'carbon dioxide',  'water', 'nitrogen' ] ),
( 'glycerol',         'calcium nitrate',    [ 'calcium oxide',       'carbon monoxide', 'water', 'nitrogen' ] ),
( 'glycerol',         'calcium nitrate',    [ 'calcium oxide',       'carbon dioxide',  'water', 'nitrogen' ] ),
( 'glycerol',         'calcium nitrate',    [ 'calcium carbonate',   'carbon dioxide',  'water', 'nitrogen' ] ),
( 'citric acid',      'calcium nitrate',    [ 'calcium oxide',       'carbon monoxide', 'water', 'nitrogen' ] ),
( 'citric acid', 	  'calcium nitrate',    [ 'calcium oxide',       'carbon dioxide',  'water', 'nitrogen' ] ),
( 'citric acid',      'calcium nitrate',    [ 'calcium carbonate',   'carbon dioxide',  'water', 'nitrogen' ] ),

( 'aluminium',        'potassium permanganate', [ 'potassium aluminate', 'aluminium oxide', 'manganese'     ] ),
( 'aluminium',        'potassium nitrate',      [ 'potassium aluminate', 'aluminium oxide', 'nitrogen'      ] ),
( 'aluminium',        'ammonium nitrate',       [ 'aluminium oxide', 'water', 'nitrogen' ] ),
( None,               'ammonium nitrate',       ['water', 'nitrogen', 'oxygen'           ] ),
( 'amonium chloride', 'ammonium nitrate',       ['water', 'nitrogen', 'hydrochloric acid'] ),
( 'carbon',           'ammonium nitrate',       [ 'carbon monoxide', 'water', 'nitrogen' ] ),
( 'carbon',           'ammonium nitrate',       [ 'carbon dioxide',  'water', 'nitrogen' ] ),
( 'sucrose',          'ammonium nitrate',       [ 'carbon monoxide', 'water', 'nitrogen' ] ),
( 'sucrose',          'ammonium nitrate',       [ 'carbon dioxide',  'water', 'nitrogen' ] ),

]


fout  = open( 'rocket_fuels.csv','w')

for reaction in reactions:
	fuel     = reaction[0]
	oxidizer = reaction[1]
	if  fuel is None:
		names = [oxidizer] + reaction[2] 
	else:
		names = [fuel] + [oxidizer] + reaction[2] 
	results = ch.reactionEnthalpy( names, echo=False )
	s = str(fuel) + '; ' + oxidizer + '; ' + ch.list2csv( results, separator='; ' )
	fout.write( s + '\n')


fout.close()
