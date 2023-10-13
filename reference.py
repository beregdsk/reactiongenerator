# rerefence heats of formation in kcal/mol from ATcT 

reference = {
	'm1_00_H2_': (0.0, 0.0, -0.0), # ATcT
	'm1_00_CH4_': (-17.812, 0.013, -0.013), # ATcT
	'm1_00_C2H6_': (-20.069, 0.031, -0.031), # ATcT
	'm1_00_C2H4_': (12.514, 0.029, -0.029), # ATcT
	'm1_00_C2H2_': (54.555, 0.031, -0.031), # ATcT 
	'm1_00_C6H6_': (19.861, 0.055, -0.055), # ATcT
	'm1_00_C2_': (197.622, 0.062, -0.062), # ATcT
	'm1_00_463-49-00_C3H4_': (45.377, 0.060, -0.060), # ATcT, allene 
	'm1_00_100-41-40_C8H10_': (7.101, 0.131, -0.131), # ATcT
	'm1_00_100-42-50_C8H8_': (35.449, 0.131, -0.131), # ATcT 
	'm1_00_106-97-80_C4H10_': (-30.062, 0.062, -0.062), # ATcT
	'm1_00_106-98-90_C4H8_': (0.012, 0.091, -0.091), # ATcT
	'm1_00_106-99-00_C4H6_': (26.489, 0.088, -0.088), # ATcT
	'm1_00_109-66-00_C5H12_': (-34.962, 0.074, -0.074), # ATcT
	'm1_00_110-82-70_C6H12_': (-29.428, 0.086, -0.086), # ATcT
	'm1_00_115-07-10_C3H6_': (4.763, 0.050, -0.050), # ATcT
	'm1_00_115-11-70_C4H8_': (-4.075, 0.100, -0.100), # ATcT
	'm1_00_142-29-00_C5H8_': (8.387, 0.110, -0.110), # ATcT ###       
	'm1_00_1828-89-32_C6H4_': (123.924, 0.287, -0.287), # ATcT     
	'm1_00_2143-69-30_C2H2_': (98.509, 0.079, -0.079), # ATcT      
	'm1_00_287-23-00_C4H8_': (6.649, 0.098, -0.098), # ATcT
	'm1_00_3355-34-80_C6H4_': (138.050, 0.287, -0.287), # ATcT
	'm1_00_4218-50-22_C2H4_': (87.612, 0.208, -0.208), # ATcT
	'm1_00_460-12-80_C4H2_': (109.955, 0.196, -0.196), # ATcT
	'm1_00_462-80-60_C6H4_': (109.730, 0.210, -0.210), # ATcT
	'm1_00_463-82-10_C5H12_': (-40.036, 0.096, -0.096), # ATcT
	'm1_00_503-17-30_C4H6_': (34.986, 0.139, -0.139), # ATcT
	'm1_00_536-74-30_C8H6_': (75.813, 0.263, -0.263), # ATcT
	'm1_00_542-92-70_C5H6_': (31.957, 0.158, -0.158), # ATcT ###
	'm1_00_594-11-60_C4H8_': (5.856, 0.311, -0.311), # ATcT
	'm1_00_624-64-60_C4H8_': (-2.672, 0.098, -0.098), # ATcT
    'm1_00_74-98-60_C3H8_': (-25.103, 0.043, -0.043), # ATcT
	'm1_00_74-99-70_C3H4_': (44.391, 0.057, -0.057), # ATcT
	'm1_00_75-19-40_C3H6_': (12.870, 0.098, -0.098), # ATcT
	'm1_00_75-28-50_C4H10_': (-32.180, 0.074, -0.074), # ATcT
	'm1_00_78-78-40_C5H12_': (-36.620, 0.103, -0.103), # ATcT
	'm1_00_78-79-50_C5H8_': (18.062, 0.160, -0.160), # ATcT
	'm1_00_822-35-50_C4H6_': (38.229, 0.222, -0.222), # ATcT
	'm1_00_108-88-30_C7H8_': (11.953, 0.081, -0.081), # ATcT
	'm1_00_206-44-0-C16H10_': (69.7, 1.0, -1.0), # NIST
	'm1_00_91-20-3_C10H8_': (35.6, 0.1, -0.1), # Our calculations based on ATcT
	'm1_00_129-00-0-C16H10_': (53.7, 0.2, -0.2), # Our calculations based on our C10H8.
}

# here all double bonds and aromatic bonds are treated in the same way
hybridization = { # (nsp, nsp2, nsp3)
	'm1_00_H2_': (0, 0, 0), 
    	'm1_00_CH4_': (0, 0, 1), # ATcT
	'm1_00_C2H6_': (0, 0, 2), # ATcT
	'm1_00_C2H4_': (0, 2, 0), # ATcT
    	'm1_00_C2H2_': (2, 0, 0), # ATcT
	'm1_00_C6H6_': (0, 6, 0), # ATcT
    	'm1_00_C2_': (0, 2, 0), # ATcT
    	'm1_00_463-49-00_C3H4_': (1, 2, 0), # ATcT, allene
	'm1_00_100-41-40_C8H10_': (0, 6, 2), # ATcT
	'm1_00_100-42-50_C8H8_': (0, 8, 0), # ATcT
	'm1_00_106-97-80_C4H10_': (0, 0, 4), # ATcT
	'm1_00_106-98-90_C4H8_': (0, 2, 2), # ATcT
	'm1_00_106-99-00_C4H6_': (0, 4, 0), # ATcT
    	'm1_00_109-66-00_C5H12_': (0, 0, 5), # ATcT
    	'm1_00_110-82-70_C6H12_': (0, 0, 6), # ATcT
	'm1_00_115-07-10_C3H6_': (0, 2, 0), # ATcT
	'm1_00_115-11-70_C4H8_': (0, 2, 2), # ATcT
	'm1_00_142-29-00_C5H8_': (0, 2, 3), # ATcT        
    	'm1_00_1828-89-32_C6H4_': (0, 6, 0), # ATcT     
    	'm1_00_2143-69-30_C2H2_': (0, 2, 0), # ATcT      
    	'm1_00_287-23-00_C4H8_': (0, 0, 4), # ATcT
    	'm1_00_3355-34-80_C6H4_': (0, 6, 0), # ATcT
    	'm1_00_4218-50-22_C2H4_': (0, 0, 2), # ATcT
    	'm1_00_460-12-80_C4H2_': (4, 0, 0), # ATcT
    	'm1_00_462-80-60_C6H4_': (2, 4, 0), # ATcT
	'm1_00_463-82-10_C5H12_': (0, 0, 5), # ATcT
    	'm1_00_503-17-30_C4H6_': (2, 0, 2), # ATcT
    	'm1_00_536-74-30_C8H6_': (2, 6, 0), # ATcT
	'm1_00_542-92-70_C5H6_': (0, 4, 1), # ATcT
    	'm1_00_594-11-60_C4H8_': (0, 0, 4), # ATcT
	'm1_00_624-64-60_C4H8_': (0, 2, 2), # ATcT
	'm1_00_74-98-60_C3H8_': (0, 0, 3), # ATcT
    	'm1_00_74-99-70_C3H4_': (2, 0, 1), # ATcT
    	'm1_00_75-19-40_C3H6_': (0, 0, 3), # ATcT
    	'm1_00_75-28-50_C4H10_': (0, 0, 4), # ATcT
    	'm1_00_78-78-40_C5H12_': (0, 0, 5), # ATcT
	'm1_00_78-79-50_C5H8_': (0, 4, 1), # ATcT
    	'm1_00_822-35-50_C4H6_': (0, 2, 2), # ATcT
    	'm1_00_108-88-30_C7H8_': (0, 6, 1), # ATcT
	'm1_00_206-44-0-C16H10_': (0,16,0), # NIST
	'm1_00_91-20-3_C10H8_': (0,10,0), # NIST
	'm1_00_129-00-0-C16H10_': (0,16,0), # NIST
}

# here all double bonds and aromatic bonds are treated in the same way
hybridization_1 = { # (nsp, nsp2, naromatic, nsp3)
    	'm1_00_H2_': (0, 0, 0, 0), 
    	'm1_00_CH4_': (0, 0, 0, 1), # ATcT
	'm1_00_C2H6_': (0, 0, 0, 2), # ATcT
	'm1_00_C2H4_': (0, 2, 0, 0), # ATcT
    	'm1_00_C2H2_': (2, 0, 0, 0), # ATcT
	'm1_00_C6H6_': (0, 0, 6, 0), # ATcT
    	'm1_00_C2_': (0, 2, 0, 0), # ATcT
    	'm1_00_463-49-00_C3H4_': (1, 2, 0, 0), # ATcT, allene
	'm1_00_100-41-40_C8H10_': (0, 0, 6, 2), # ATcT
	'm1_00_100-42-50_C8H8_': (0, 2, 6, 0), # ATcT
	'm1_00_106-97-80_C4H10_': (0, 0, 0, 4), # ATcT
	'm1_00_106-98-90_C4H8_': (0, 2, 0, 2), # ATcT
	'm1_00_106-99-00_C4H6_': (0, 4, 0, 0), # ATcT
    	'm1_00_109-66-00_C5H12_': (0, 0, 0, 5), # ATcT
    	'm1_00_110-82-70_C6H12_': (0, 0, 0, 6), # ATcT
	'm1_00_115-07-10_C3H6_': (0, 2, 0, 1), # ATcT
	'm1_00_115-11-70_C4H8_': (0, 2, 0, 2), # ATcT
	'm1_00_142-29-00_C5H8_': (0, 2, 0, 3), # ATcT        
    	'm1_00_1828-89-32_C6H4_': (0, 0, 6, 0), # ATcT     
    	'm1_00_2143-69-30_C2H2_': (0, 2, 0, 0), # ATcT      
    	'm1_00_287-23-00_C4H8_': (0, 0, 0, 4), # ATcT
    	'm1_00_3355-34-80_C6H4_': (0, 0, 6, 0), # ATcT
    	'm1_00_4218-50-22_C2H4_': (0, 0, 0, 2), # ATcT
    	'm1_00_460-12-80_C4H2_': (4, 0, 0, 0), # ATcT
    	'm1_00_462-80-60_C6H4_': (2, 4, 0, 0), # ATcT
	'm1_00_463-82-10_C5H12_': (0, 0, 0, 5), # ATcT
    	'm1_00_503-17-30_C4H6_': (2, 0, 0, 2), # ATcT
    	'm1_00_536-74-30_C8H6_': (2, 0, 6, 0), # ATcT
	'm1_00_542-92-70_C5H6_': (0, 4, 0, 1), # ATcT
    	'm1_00_594-11-60_C4H8_': (0, 0, 0, 4), # ATcT
	'm1_00_624-64-60_C4H8_': (0, 2, 0, 2), # ATcT
	'm1_00_74-98-60_C3H8_': (0, 0, 0, 3), # ATcT
    	'm1_00_74-99-70_C3H4_': (2, 0, 0, 1), # ATcT
    	'm1_00_75-19-40_C3H6_': (0, 0, 0, 3), # ATcT
    	'm1_00_75-28-50_C4H10_': (0, 0, 0, 4), # ATcT
    	'm1_00_78-78-40_C5H12_': (0, 0, 0, 5), # ATcT
	'm1_00_78-79-50_C5H8_': (0, 4, 0, 1), # ATcT
    	'm1_00_822-35-50_C4H6_': (0, 2, 0, 2), # ATcT
    	'm1_00_108-88-30_C7H8_': (0, 0, 6, 1), # ATcT
	'm1_00_206-44-0-C16H10_': (0, 0, 16, 0), # NIST
	'm1_00_91-20-3_C10H8_': (0,0,10,0), # NIST
	'm1_00_129-00-0-C16H10_': (0,0,16,0), # NIST
}

name = {
	'm1_00_H2_': 'Hydrogen', # ATcT
	'm1_00_CH4_': 'Methane', # ATcT
	'm1_00_C2H6_': 'Ethane', # ATcT
	'm1_00_C2H4_': 'Ethylene', # ATcT
	'm1_00_C2H2_': 'Acetylene', # ATcT
	'm1_00_C6H6_': 'Benzene', # ATcT
	'm1_00_C2_': 'Ethynylene', # ATcT
	'm1_00_463-49-00_C3H4_': 'Allene', # ATcT, allene
	'm1_00_100-41-40_C8H10_': 'Ethylbenzene', # ATcT
	'm1_00_100-42-50_C8H8_': 'Phenylethene', # ATcT
	'm1_00_106-97-80_C4H10_': 'n-Butane', # ATcT
	'm1_00_106-98-90_C4H8_': '1-Butene', # ATcT
	'm1_00_106-99-00_C4H6_': '1,3-Butadiene', # ATcT
	'm1_00_109-66-00_C5H12_': 'n-Pentane', # ATcT
	'm1_00_110-82-70_C6H12_': 'Cyclohexane', # ATcT
	'm1_00_115-07-10_C3H6_': 'Propene', # ATcT
	'm1_00_115-11-70_C4H8_': 'Isobutene', # ATcT
	'm1_00_142-29-00_C5H8_': 'Cyclopentene', # ATcT        
	'm1_00_1828-89-32_C6H4_': 'm-Benzyne', # ATcT     
	'm1_00_2143-69-30_C2H2_': 'Vinylidene', # ATcT      
	'm1_00_287-23-00_C4H8_': 'Cyclobutane', # ATcT
	'm1_00_3355-34-80_C6H4_': 'p-Benzyne', # ATcT
	'm1_00_4218-50-22_C2H4_': 'Ethylidene', # ATcT
	'm1_00_460-12-80_C4H2_': '1,3-Butadiyne', # ATcT
	'm1_00_462-80-60_C6H4_': 'o-Benzyne', # ATcT
	'm1_00_463-82-10_C5H12_': 'neo-Pentane', # ATcT
	'm1_00_503-17-30_C4H6_': '2-Butyne', # ATcT
	'm1_00_536-74-30_C8H6_': 'Phenylacetylene', # ATcT
	'm1_00_542-92-70_C5H6_': 'Cyclopentadiene', # ATcT
        'm1_00_594-11-60_C4H8_': 'Methylcyclopropane', # ATcT
	'm1_00_624-64-60_C4H8_': 'trans-2-Butene', # ATcT
        'm1_00_74-98-60_C3H8_': 'Propane', # ATcT
        'm1_00_74-99-70_C3H4_': 'Propyne', # ATcT
	'm1_00_75-19-40_C3H6_': 'Cyclopropane', # ATcT
	'm1_00_75-28-50_C4H10_': 'iso-Butane', # ATcT
	'm1_00_78-78-40_C5H12_': 'iso-Pentane', # ATcT
	'm1_00_78-79-50_C5H8_': 'Isoprene', # ATcT
	'm1_00_822-35-50_C4H6_': 'Cyclobutene', # ATcT
	'm1_00_108-88-30_C7H8_': 'Toluene', # ATcT
	'm1_00_206-44-0-C16H10_': 'Fluoranthene', # NIST
	'm1_00_91-20-3_C10H8_': 'Naphthalene', # NIST
	'm1_00_129-00-0-C16H10_': 'Pyrene', # NIST
}