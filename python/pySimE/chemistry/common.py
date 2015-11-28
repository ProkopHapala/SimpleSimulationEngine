#!/usr/bin/python

import re
import numpy as np

from .data.substances import substance_dict, substance_param_names
from .data.elements   import element_dict, element_param_names

# ======== parse reaction

def str2composition( s ):
	d = {}
	for e,m in re.findall('([A-Z][a-z]?)([0-9]*)',s):
		m=1 if m=='' else int(m)
		d[e]=int(m)
	return d

def parseReactionSide( s ):
	names = []
	comps = []
	for sub in s.split('+'):
		names.append( sub )
		comps.append( str2composition( sub ) )
	return names, comps

def parseReaction( s ):
	sides = s.split('=')
	return ( ( parseReactionSide( sides[0] ) ), ( parseReactionSide( sides[1] ) ) )

def i2s_but1( i ):
	s = ''
	if i <> 1:
		s = str(i)
	return s

def composition2string( d ):
	s = ''
	for key, value in d.iteritems():
		s = s + key + i2s_but1( value )
	return s

def reaction2string( names, coefs ):
	sleft  = ''
	sright = ''
	for i,c in enumerate(coefs):
		cstr = ' + '+i2s_but1( abs(c) )+' '
		if c > 0:
			sright = sright + cstr + names[i]
		else:
			sleft  = sleft + cstr + names[i]
	return sleft + ' = ' + sright

def printNiceMap( keys, values ):
	for i,key in enumerate(keys):
		print key,' : ', values[i]

# ======== database I/O

def takeSubstance( keys, param, substance_dict=substance_dict ):
	return [ substance_dict[key] for key in keys ]

def takeSubstanceParam( keys, param, substance_dict=substance_dict, substance_param_names=substance_param_names ):
	index = substance_param_names[param]
	return [ substance_dict[key][index] for key in keys ]

def names2comps( names, substance_dict=substance_dict ):
	comps = []
	for name in names:
		comps.append( str2composition( substance_dict[name][1] ) )
	return comps

def molarMass( comp, element_dict=element_dict ):
	M = 0
	for key,value in comp.iteritems():
		M += element_dict[key][1] * value
	return M

def comps2mass( comps, element_dict=element_dict ):
	Ms = []
	for comp in comps:
		Ms.append( molarMass( comp, element_dict ) )
	return Ms

def selectPhase( names, pahse='g', substance_dict=substance_dict ):
	return [ (substance_dict[name][0]==pahse) for name in names ]


def list2csv( l, separator=';' ):
	s = ''
	for item in l:
		s = s + str(item) + ';'
	return s

# ======== atomic balance in reaction

def listElements( comps ):
	elems=set()
	for comp in comps:
		for elem in comp.keys():
			elems.add(elem)
	return list(elems)
	
def reactionMatrix( comps, elems, M=None ):
	if M==None:
		M = np.zeros( (len(elems),len(comps)) )
	for i,comp in enumerate(comps):
		for j,elem in enumerate(elems):
			if elem in comp:
				M[j,i] = comp[elem]
	return M

def atomicBalance( names, comps, pivot=0 ):
	elems   = listElements( comps )
	ncomp   = len(comps)
	nelem   = len(elems)
	nrow    = nelem
	if nelem < ncomp:
		nrow = nelem + 1
	M       = np.zeros( (nrow,ncomp) )
	M       = reactionMatrix( comps, elems, M=M );
	rhs     = np.zeros  ( len( M[0] ), dtype=np.int );  
	if nelem < ncomp:
		M[-1,pivot] = 1
		rhs[-1] = -1;
	else:
		rhs[0]  = -1;
	#print M
	#print rhs
	coefs    = np.linalg.solve( M, rhs )
	#print M
	#print "rhs: ",rhs
	#print "coefs: ",coefs
	#print np.dot(M,coefs)
	return coefs

def balanceReactionString( s, pivot=0 ):
	rl = parseReaction( s )
	names = rl[0][0] + rl[1][0]
	comps = rl[0][1] + rl[1][1]
	coefs = atomicBalance( names, comps, pivot=pivot )
	return reaction2string( names, coefs )

def reactionEnthalpy( names,   element_dict=element_dict, substance_dict=substance_dict, echo=True ):
	comps     = names2comps  ( names, substance_dict ) 
	masses    = comps2mass   ( comps, element_dict   )
	coefs     = atomicBalance( names, comps, pivot=0    )
	#coefs_left  = -np.fmin( coefs, 0 )
	coefs_right =  np.fmax( coefs, 0 )
	Hs        = takeSubstanceParam( names, 'enthalpy', substance_dict=substance_dict )
	formulas  = takeSubstanceParam( names, 'formula',  substance_dict=substance_dict )
	gasMask   = selectPhase       ( names, pahse='g',  substance_dict=substance_dict )
	Mtot      = sum( coefs_right * np.array( masses  )   )
	Htot      = sum( coefs       * np.array( Hs      )   )
	gasMoles  = sum( coefs_right[  np.array( gasMask ) ] )
	if echo:
		print reaction2string( formulas, coefs )
		print " enthalpy_change      = ", Htot,                 " kJ/mol "
		print " total mass           = ", Mtot,                 " g/mol  "
		print " Energy density       = ", (Htot / Mtot) ,       " MJ/kg  "
		print " gasMoles             = ", gasMoles,             " mol    "
		print " effective molar mass = ", ( Mtot / gasMoles ) , " g/mol  "
		print " energy gas mol       = ", ( Htot / gasMoles ) , " kJ/mol "
	return [ Htot, Mtot, gasMoles, (Htot / Mtot), ( Mtot / gasMoles ), ( Htot / gasMoles ), reaction2string( formulas, coefs ), coefs ]  




