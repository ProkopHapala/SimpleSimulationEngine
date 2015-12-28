#!/usr/bin/python

'''
REFERENCES:
http://www.cementkilns.co.uk/data_enthalpy_formation.html
http://bilbo.chm.uri.edu/CHM112/tables/thermtable.htm
http://webbook.nist.gov/chemistry/
'''


substance_param_names ={
'pahse'    : 0,
'formula'  : 1,
'enthalpy' : 2,
}


substance_dict = {
'hydrogen'		        : ['g','H2',  0 ],
'carbon'		        : ['s','C',   0 ],
'nitrogen'		        : ['g','N2',  0 ],
'oxygen'		        : ['g','O2',  0 ],
'sulphur'		        : ['s','S',   0 ],
'chlorine'              : ['g','Cl2', 0 ],
'potassium'             : ['s','K',   0 ],
'aluminium'             : ['s','Al',  0 ],
'manganese'             : ['s','Mn',  0 ],
'iron'                  : ['s','Fe',  0 ],

'methan'		        : ['g','CH4',     -74.9 ],
'acetylene'		        : ['g','C2H2',   +226.8 ],
'n-hexane'		        : ['g','C6H14',  -167.4 ],
'n-decane'		        : ['l','C10H22', -249.4 ],

'glucose'               : ['s','C6H12O6',   -1273.3 ],
'sucrose'               : ['s','C12H22O11', -2226.1 ],
'ceulose'				: ['s','C6H10O5', 	-963.0  ],
'glycerol'				: ['l','C3H8O3',    -669.6  ],
'glycerol trinitrate'   : ['l','C3H5N3O9',  -370    ],
'citric acid'			: ['s','C6H8O7', 	-1548.8 ],

'amonia'				: ['g','NH3',		-46		],
'hydrazine'				: ['g','N2H4',		+50.63	],
'nitric acid'           : ['l','HNO3',      -207    ],
'nitrous oxide'         : ['g','N2O',       +82.05  ],
'nitric oxide'          : ['g','NO',		+91.29  ],
'nitrogen dioxide'      : ['g','NO2',       +34.0   ],
'dinitrogen tetroxide'  : ['g','N2O4',      +9.16   ],
'dinitrogen pentoxide'  : ['s','N2O5',		-43.1	],

'ammonium nitrate'		: ['s','N2H4O3',  	-365.6  ],
'ammonium perchlorate'	: ['s','NH4ClO4',	-295.77	],
'amonium chloride'		: ['s','NH4Cl'	,	-314.4  ],

'hydrochloric acid'		: ['g','HCl',		-92.31	],

'water'                 : ['g','H2O',       -285.8  ],
'hydrogen peroxide'		: ['g','H2O2',		-187.80	],
'ozone'					: ['g','O3',		+142.67 ],

'carbon monoxide'       : ['g','CO',        -110.5 ],
'carbon dioxide'        : ['g','CO2',       -393.5 ],

'aluminium oxide'       : ['s','Al2O3',      -1675.7 ],
'manganese dioxide'     : ['s','MnO2',       -520    ],
'iron(III) oxide'       : ['s','Fe2O3',      -824.2  ],

'potassium nitrate'     : ['s','KNO3',    -494     ],
'potassium oxide'       : ['s','K2O',     -363     ],
'potassium hydroxide'   : ['s','KOH',     -425.9   ],
'potassium carbonate'   : ['s','K2CO3',   -1150    ],
'potassium chloride'    : ['s','KCl',     -436.68  ],
'potassium perchlorate' : ['s','KClO4',   -433     ],
'potassium permanganate': ['s','KMnO4',   -813.4   ],
'potassium aluminate'   : ['s','K2Al2O4', -2281.10 ],

'calcium nitrate'       : ['s','CaN2O6',   -938.4 ],
'calcium oxide'         : ['s','CaO',      -635.1 ],
'calcium hydroxide'     : ['s','CaO2H2',   -986.1 ],
'calcium carbonate'     : ['s','CaCO3',    -1207  ],
'calcium chloride'      : ['s','CaCl2',    -795.8  ],





}

