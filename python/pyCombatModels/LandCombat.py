
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os
import sys
#from path import path

work_dir    = os.path.dirname( os.path.realpath( __file__ ) )
python_path = os.path.normpath( work_dir  + '../../' );    print "python_path :  ", python_path
sys.path.append( python_path   );                          print "sys.path    :  ", sys.path

from   pyMeta import cpp_utils
from   pyMeta.cpp_utils import _np_as
#import cpp_utils 
#from cpp_utils import _np_as

c_double_p = ctypes.POINTER(c_double)

# ===== To generate Interfaces automatically from headers call:
header_strings = [
    "int newWeaponType( char* str ){",
    "int newUnitType( char* str ){",
    "int newUnit( int ityp ){",
    "void setUnitPose( int i, double* pos, double* rot ){",
    "double shoot( int itarget, int iattacker, int igun ){"
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_name='LandCombatLib'
cpp_utils.make(cpp_name)
lib  = cpp_utils.loadLib( cpp_name )

# =============== C / Python interfaces 

#  int newWeaponType( const char* str ){
lib.newWeaponType.argtypes  = [c_char_p] 
lib.newWeaponType.restype   =  c_int
def newWeaponType(str):
    return lib.newWeaponType( _np_as(str,c_char_p) ) 

#  int newUnitType( const char* str ){
lib.newUnitType.argtypes  = [c_char_p] 
lib.newUnitType.restype   =  c_int
def newUnitType(str):
    return lib.newUnitType( _np_as(str,c_char_p) ) 

#  int newUnit( int ityp ){
lib.newUnit.argtypes  = [c_int] 
lib.newUnit.restype   =  c_int
def newUnit(ityp):
    return lib.newUnit(ityp) 

#  void setUnitPose( int i, double* pos, double* rot ){
lib.setUnitPose.argtypes  = [c_int, c_double_p, c_double_p] 
lib.setUnitPose.restype   =  None
def setUnitPose(i, pos, rot):
    return lib.setUnitPose(i, pos.ctypes.data_as(c_double_p), rot.ctypes.data_as(c_double_p)) 

#  double shoot( int itarget, int iattacker, int igun ){
lib.shoot.argtypes  = [c_int, c_int, c_int] 
lib.shoot.restype   =  c_double
def shoot(itarget, iattacker, igun):
    return lib.shoot(itarget, iattacker, igun) 

# =============== Test Run

if __name__ == "__main__":  
    import matplotlib.pyplot as plt

    # see
    # /cpp/apps/LandTactics/data/GunTypes.ini
    # /cpp/apps/LandTactics/data/UnitTypes.ini

    newWeaponType( "granade	;		0.03	1	0.02		10	0.005	;	0.003	2.00E+05	;	Frag" )
    newWeaponType( "SMG9mm	;		0.2	10	1.00E-02	300	0.01	;	0.005	0.00E+00	;	KE" )
    newWeaponType( "RG7.62mm	;		0.1	1	1.00E-03	1000	0.01	;	0.015	0.00E+00	;	KE" )
    newWeaponType( "MG7.62mm;		0.1	10	5.00E-03	1000	0.01	;	0.015	0.00E+00	;	KE" )
    newWeaponType( "MG12.7mm;		0.1	10	5.00E-03	1000	0.043	;	0.03	0.00E+00	;	KE" )
    newWeaponType( "KwK2cm	;		0.1	10	5.00E-03	1000	0.12	;	0.04	0.00E+00	;	KE" )
    newWeaponType( "PanzerFaust;		0.1	1	1.00E-03	60	7	;	0.25	0.00E+00	;	HEAT" )
    newWeaponType( "AT37mm	;		0.1	1	1.00E-03	1000	0.5	;	0.05	1.00E+04	;	KE" )
    newWeaponType( "AT50mm	;		0.1	1	1.00E-03	1000	2	;	0.07	1.00E+05	;	KE" )
    newWeaponType( "AT75mm	;		0.1	1	1.00E-03	1000	4	;	0.15	2.00E+05	;	KE" )
    newWeaponType( "AT88mm	;		0.1	1	1.00E-03	1000	8	;	0.2	2.00E+05	;	KE" )

    infSMG    = newUnitType  ( " SMG		;	inf	;	1	1	2	;	80		1.5	1.00E+03;	0	0	0	0	0	;	1	;	2	;	SMG9mm ;	granade ") 
    inf       = newUnitType  ( " Rifle		;	inf	;	1	1	2	;	80		1.3	1.00E+03;	0	0	0	0	0	;	1	;	2	;	RG7.62mm ;	granade ") 
    infMG     = newUnitType  ( " MG		;	inf	;	1	1	2	;	80		1.1	1.00E+03;	0	0	0	0	0	;	1	;	1	;	MG7.62mm	") 
    infAT     = newUnitType  ( " PanzerFaust	;	inf	;	1	1	2	;	80		1.1	1.00E+03;	0	0	0	0	0	;	1	;	1	;	PanzerFaust ") 	
    atg37     = newUnitType  ( " PAK_37mm	;	gun	;	3	2	1.2	;	370		0.7	1.00E+03;	5	0	0	0	0	;	3	;	1	;	AT37mm	") 
    atg50     = newUnitType  ( " PAK_50mm	;	gun	;	4	2	1.4	;	830		0.5	1.00E+03;	5	0	0	0	0	;	4	;	1	;	AT50mm	") 
    atg       = rifle = newUnitType  ( " PAK_75mm	;	gun	;	5	2	1.6	;	1400		0.3	1.00E+03;	5	0	0	0	0	;	5	;	1	;	AT75mm	") 
    tPz2      = newUnitType  ( " PzII_tank	;	tank	;	5	2.2	2	;	8900		15	1.00E+03;	35	15	15	20	5	;	4	;	2	;	KwK2cm ;	MG7.62mm") 
    tM4       = newUnitType  ( " M4A1_tank	;	tank	;	6	3	3	;	3.50E+04	12	1.00E+03;	80	40	30	20	20	;	4	;	2	;	AT75mm ;	MG12.7mm") 
    tPz4rifle = newUnitType  ( " PzIV_tank	;	tank	;	6	2.9	2.6	;	2.50E+04	13	1.00E+03;	80	30	30	20	20	;	4	;	2	;	AT75mm ;	MG7.62mm") 
    tPanther  = newUnitType  ( " Panther_tank	;	tank	;	6.8	3.5	3	;	5.40E+04	10	1.00E+03;	150	40	16	20	20	;	4	;	2	;	AT88mm ;	MG7.62mm") 

    u1 = newUnit( infAT )
    u2 = newUnit( tM4   )

    dmg = shoot( u1, u2, 1 )
    print( "dmg ", dmg )

    #newWeaponType(str)
    #newUnitType  (str)

    #plt.show()