#!/usr/bin/python

import pyMeta as meta

cpp_name = "/home/prokop/git/SimpleSimulationEngine/cpp/libs/Molecular/Molecular.cpp"
py_name  = "/home/prokop/git/SimpleSimulationEngine/python/c_auto_interface_Molecular.py"

func_names = [
"getPlaneWaveDescriptor",
"compDistanceT",
"initComparatorT",
"compDistance",
"initComparator",
"testMultipole",
"collisionForce"
]


#with open(cpp_name) as f:
#    cpp_lines = f.readlines()
#    meta.findCloseBrace( cpp_lines, 0, 0, q='{', p='}' )

print meta.prim_types.values()
meta.generate_interface( cpp_name, py_name, func_names )

