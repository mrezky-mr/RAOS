#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 17 20:41 2022
Updated on -

@author: M. Rezky
@PI: Endang Soegiartini

input.py was developed for generating all 243 (3^5) clones of asteroid orbit.

This program are needed as sim.add() command that links with HORIZONS API
apparently called the ephimerides calculation 
(https://ssd.jpl.nasa.gov/horizons/app.html#/), not the values from observation 
results (https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/). Please copy all 
numbers from the latter link, as it was based on observation logs.

"""
########### OBJECT INFORMATION INPUT STARTS HERE #################

obj_name = '2022BJ8'
date = '2022-08-09 00:00'

'''Semi-major axis distance (asmd), a [au]'''
asmd_value = 0.7852852601007706
asmd_sigma = 1.327E-6

'''Orbital eccentricity (eccs), e [none]'''
eccs_value = 0.2486899450820991
eccs_sigma = 6.4422E-6

'''Orbital inclination (incl), i [deg]'''
incl_value = 15.82521869196891
incl_sigma = .00034608

'''Longitude of ascending node (ascn), OMEGA [deg]'''
ascn_value = 95.84694910204503
ascn_sigma = .00032108

'''Argument of perihelion (peri), omega [deg]'''
peri_value = 293.872626469623
peri_sigma = .00024355


########### OBJECT INFORMATION INPUT ENDS HERE ###################
##################################################################


import numpy as np
import pandas as pd
from datetime import datetime
from itertools import product as prd

st_time = datetime.now()
start = st_time.strftime("%Y%m%d_%H%M%S")
print(start)

asmd_v = [asmd_value, asmd_value+asmd_sigma, asmd_value-asmd_sigma]
eccs_v = [eccs_value, eccs_value+eccs_sigma, eccs_value-eccs_sigma]
incl_v = np.radians([incl_value, incl_value+incl_sigma, incl_value-incl_sigma])
ascn_v = np.radians([ascn_value, ascn_value+ascn_sigma, ascn_value-ascn_sigma])
peri_v = np.radians([peri_value, peri_value+peri_sigma, peri_value-peri_sigma])

'''
Input parameter set randomizers
[0] = nominal value
[1] = nominal + 1-sigma value
[2] = nominal - 1-sigma value
'''
clist = list(prd([0,1,2], repeat=5))

numb = []
asmd = []
eccs = []
incl = []
ascn = []
peri = []
for i in range(len(clist)):
    numb.append(i)
    asmd.append(asmd_v[clist[i][0]])
    eccs.append(eccs_v[clist[i][1]])
    incl.append(incl_v[clist[i][2]])
    ascn.append(ascn_v[clist[i][3]])
    peri.append(peri_v[clist[i][4]])
    
dict = {"index": numb, "asmd": asmd, "eccs": eccs, "incl": incl, 
        "ascn": ascn, "peri": peri}
df = pd.DataFrame(dict)
df.to_csv('params_' + obj_name + '.csv')

file = open("input_" + obj_name + "_" + start + ".txt","w")
file.write("SIMULATION INPUT SUMMARY \n")
file.write("\n")
file.write("Object name                     : {a} \n".format(a=obj_name))
file.write("Ephimerides Epoch on HORIZONS   : {b} \n".format(b=date))
file.write("---------- \n")
file.write("\n")
file.write("ORBITAL PARAMETERS, NOMINAL VALUE \n")
file.write("Semi-major axis distance (au)   : {a} \n".format(a=asmd_value))
file.write("Eccentricity                    : {b} \n".format(b=eccs_value))
file.write("Inclination (deg)               : {a} \n".format(a=incl_value))
file.write("Longitude, ascending node (deg) : {a} \n".format(a=ascn_value))
file.write("Argument of perihelion (deg)    : {a} \n".format(a=peri_value))
file.write("\n")
file.write("ORBITAL PARAMETERS, 1-SIGMA VALUE \n")
file.write("Semi-major axis distance (au)   : {a} \n".format(a=asmd_sigma))
file.write("Eccentricity                    : {b} \n".format(b=eccs_sigma))
file.write("Inclination (deg)               : {a} \n".format(a=incl_sigma))
file.write("Longitude, ascending node (deg) : {a} \n".format(a=ascn_sigma))
file.write("Argument of perihelion (deg)    : {a} \n".format(a=peri_sigma))
file.write("---------- \n")
file.close()