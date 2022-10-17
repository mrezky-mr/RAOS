#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 11 2022 
Update on Oct 17 2022

@author: M. Rezky
@PI: E. Soegiartini

N-body simulation program for all 243 asteroid orbit clones, with only gravity 
forces involved. The gravity calculation was based on the asteroid's 
interactions with Sun, 8 planets, and the Moon.

Please review all input parameters before running this program. Make sure that
the input parameter file ('params_{asteroidname}.csv') and output files are
located in the same directory (the one included in 'datadir' input parameter
below).
"""

from __future__ import print_function

# SIMULATION INPUT ###########################################################
# INTEGRATOR SETTINGS
integrator = "ias15"
grav_mode = "compensated"
int_tolerance = 1e-9
gr_type = "gr_potential"             

# COLLISION AND BOUNDARY OPTIONS
col_detmode = "none"
col_resmode = "halt"
sim_time = 1e6              # simulation time limit [yr]
sim_dast = 10               # simulation data timestep [yr]

# TARGET OBJECT INFORMATION
obj = '2021PH27'
date = '2022-08-09 00:00'
datadir = '~/atr/atr2/'

##############################################################################

'''REQUIRED MODULES'''
import rebound
import reboundx
import math as m
import numpy as np
import pandas as pd
from datetime import datetime

'''TIMESTAMP OF THE START OF THE PROGRAM'''
st_time = datetime.now()
start = st_time.strftime("%Y%m%d_%H%M%S")
print("sim_grpn.py starts running at ", start)
print('-----')

'''CONVERSION VARIABLES'''
au = 1.495978707e11             # astronomical units [m]
yr = 365.25*24*3600             # Julian year [s]
jc = 100                        # Julian century [yr]
arcsec = 4.84813681109536e-06   # arcsec [rad]
lightsp = 299792458 * yr/au     # speed of light [au/yr]

'''HEARTBEAT FUNCTION FOR NUMERICAL INTEGRATION & DATA COLLECTION'''
def heartbeat(sim):
    global ER0, EH0, AM0, n, simx
    smr = sim.contents
    sp2 = smr.particles
    t = smr.t
    qt = 0
    
    if t > n * sim_dast:
        for q in range(smr.N - smr.N_active):
            labl = 'A{:03d}'.format(q)
            
            '''
            CHECK IF THIS CLONE HAVE COLLIDED WITH SUN/PLANET/MOON BEFORE
            IF YES, THIS CLONE IS SKIPPED FOR THE REST OF THE SIMULATION
            '''
            saved_col = False
            for u in c:
                if u == labl:
                    saved_col = True
                    break
                else:
                    saved_col = False
            
            if saved_col == False:
                '''INITIATE ARRAYS FOR DATA COLLECTION'''
                ## ORBITAL PARAMETERS OF TARGET OBJECT
                tis = []    # simulation time [yr]
                ecc = []    # eccentricity (e)
                asm = []    # semi-major distance (a) [au]
                inc = []    # inclination (i) [deg]
                lan = []    # longitude of ascending node (Omega) [deg]
                aop = []    # argument of perihelion (omega) [deg]
                ftr = []    # true anomaly (f) [deg]
                ## SYSTEM DIAGNOSTIC PARAMETERS
                erg = []    # REBOUND energy
                err = []    # REBOUND relative energy error
                anm = []    # REBOUND angular momentum
                anr = []    # REBOUND relative am error
                ehg = []    # Hamiltonian energy
                ehr = []    # Hamiltonian relative energy error
                ## RESONANCE ANALYSIS PARAMETERS (W/ MERCURY & VENUS)
                amm = []    # semi-major distance (a) of Mercury [au]
                ecm = []    # eccentricity (e) of Mercury
                aom = []    # argument of perihelion (omega) of Mercury [deg]
                amv = []    # semi-major distance (a) of Venus [au]
                ecv = []    # eccentricity (e) of Venus
                aov = []    # semi-major distance (a) of Venus [deg]
                ## GENERAL RELATIVITY DIAGNOSTIC PARAMETERS
                pmm = []    # Mercury pericenter longitude (pomega) [rad]
                pmc = []    # Mercury rate of perihelion precession ["/JC]
                pla = []    # target object's pericenter longitude (pomega) [rad]
                plc = []    # target object's rate of perihelion precession ["/JC]
                
                tis.append(t)
                
                '''CALLING PERIHELION PRECISION DATA FROM LAST TIME STEP'''
                if n != 0:
                    datacol = pd.read_csv(datadir + 'grpn_' + start + '_data_' + labl + '.csv')
                    tisl = datacol['t']
                    pmml = datacol['gre_ma']
                    plal = datacol['gre_oa']
                    deltat = t - tisl[n-1]
    
                '''ENERGY CALCULATION'''
                enr = smr.calculate_energy()
                enh = simx.gr_hamiltonian(gr)
                if n == 0:
                    ER0 = enr
                    EH0 = enh
                    rer = 0
                    reh = 0
                else:
                    rer = (enr - ER0) / ER0
                    reh = (enh - EH0) / EH0
                erg.append(enr)
                err.append(rer)
                ehg.append(enh)
                ehr.append(reh)
    
                '''ANGULAR MOMENTUM CALCULATION'''
                amx, amy, amz = smr.calculate_angular_momentum()
                ama = m.sqrt(amx**2 + amy**2 + amz**2)
                if n == 0:
                    AM0 = ama
                    amd = 0
                else:
                    amd = (ama - AM0) / AM0
                anm.append(ama)
                anr.append(amd)
        
                '''MERCURY PERIHELION PRECESSION CALCULATION (for control)'''  
                pm = sp2[1].pomega
                if n == 0:
                    deltpm = 0
                    pmr = 0
                else:                   
                    pml = pmml[n-1]
                    deltpm = pm - pml
                    pmr = deltpm * jc / deltat / arcsec
                pmm.append(pm)
                pmc.append(pmr)
                    
                '''CLOSE ENCOUNTER DETECTION'''
                for body_i in range(len(main_bodies)):
                    if body_i != 0:
                        ce_rildis = sp2[body_i] ** sp2[labl]
                        ce_mindis = sp2[body_i].rhill   # simplified formula, check "tools.c"
                        if ce_rildis <= ce_mindis:
                            dice = {"Asteroid clone": [labl], "Object name": [pname[body_i]], 
                                    "CE_time": [t], "Distance [au]": [ce_rildis], 
                                    "Distance [R_Hill]": [ce_rildis/ce_mindis]}
                        else:
                            dice = {}
                        df = pd.DataFrame(dice)
                        if n == 0:
                            df.to_csv('grpn_' + start + '_ce-data.csv')
                        else:
                            df.to_csv('grpn_' + start + '_ce-data.csv', mode='a', header=False)
                
                '''ASTEROID PERIHELION PRECESSION CALCULATION'''  
                po = sp2[labl].pomega
                if n == 0:
                    deltpo = 0
                    por = 0
                else:
                    pol = plal[n-1]               
                    deltpo = po - pol
                    por = deltpo * jc / deltat / arcsec
                pla.append(po)
                plc.append(por)
        
                '''ORBITAL ELEMENT DATA COLLECTION'''
                ecc.append(sp2[labl].e)                                  
                asm.append(sp2[labl].a)                                  
                inc.append(m.degrees(sp2[labl].inc))                   
                lan.append(m.degrees(sp2[labl].Omega))                   
                aop.append(m.degrees(sp2[labl].omega))                        
                ftr.append(m.degrees(sp2[labl].f))
        
                '''RESONANCE ANALYSIS DATA COLLECTION'''
                amm.append(sp2[1].a)
                ecm.append(sp2[1].e)
                aom.append(m.degrees(sp2[1].omega))
                amv.append(sp2[2].a)   
                ecv.append(sp2[2].e)   
                aov.append(m.degrees(sp2[2].omega))
                
                '''CHECKING COLLISION BETWEEN THE ASTEROID CLONE AND EVERY ACTIVE OBJECT'''
                for k in range(smr.N_active):
                    if k != 0:
                        rr = sp2[k].r
                        dx = abs(sp2[k].x - sp2[labl].x)
                        dy = abs(sp2[k].y - sp2[labl].y)
                        dz = abs(sp2[k].z - sp2[labl].z)
                        dr = np.sqrt(dx**2 + dy**2 + dz**2)
                    else:
                        rr = 0.00465047
                        dr = sp2[labl].d
                    if dr <= rr:
                        file.write('Clone ' + labl + ' collided with ' + pname[k] + ' on t=' + str(round(t,2)))
                        c.append(labl)
            
                '''STORING ALL CALCULATED PARAMETERS AS A CSV'''
                dict = {"t": tis, "obj_e": ecc, "obj_a": asm, "obj_i": inc, 
                        "obj_O": lan, "obj_o": aop, "obj_f": ftr, "erg_ra": erg, 
                        "erg_re": err, "erg_ha": ehg, "erg_he": ehr, "anm_a": anm, 
                        "anm_e": anr, "res_ma": amm, "res_me": ecm, "res_mo": aom, 
                        "res_va": amv, "res_ve": ecv, "res_vo": aov, "gre_oa": pla, 
                        "gre_oc": plc, "gre_ma": pmm, "gre_mc": pmc}
                df = pd.DataFrame(dict)
                if n == 0:
                    df.to_csv('grpn_' + start + '_data_' + labl + '.csv')
                else:
                    df.to_csv('grpn_' + start + '_data_' + labl + '.csv', mode='a', header=False)
            
            '''LIVE UPDATE OF THE SIMULATION STATUS'''
            if q == qt and saved_col == False:
                print('-----')
                if qt == 0:
                    print('Results from nominal asteroid orbit')
                else:
                    print('Results from clone ' + str(qt) + ' asteroid orbit')
                print("Integration time (yr)        : ", t)
                print("Rel. energy error            : ", rer)
                print("Rel. Hamiltonian error       : ", reh)
                print("Rel. angular momentum error  : ", amd)
                print('Obj. pomega chg. rate ("/JC) : ', por)
                print('Mer. pomega chg. rate ("/JC) : ', pmr)
                print("Clones collided              : ", len(c))
            elif q == qt and saved_col == True:
                qt += 1
        n += 1
    
'''INITIATE REBOUND AND THE DEFAULT UNITS'''
sim = rebound.Simulation()
sim.units = ('yr','AU','kg')
ps = sim.particles

'''SUBMIT INTEGRATOR & COLLISION SETTINGS'''
sim.integrator = integrator
sim.gravity = grav_mode
if integrator == 'ias15':
    sim.ri_ias15.epsilon = int_tolerance
sim.collision = col_detmode
sim.collision_resolve = col_resmode
sim.dt = 1/365.25

'''ADDING PLANETS AND TARGET OBJECT (DEFAULT: SUN, 8 PLANETS and MOON/301)'''
main_bodies = ["10","199","299","399","301","499","599","699","799","899"]
sim.add(main_bodies, date=date)
sim.N_active = sim.N
sim.testparticle_type = 0       # put '1' if there asteroid mass is available.

## TEST PARTICLE STATUS ARE NOT CARRIED AWAY BY REBOUNDX
## CLONES ARE SUBJECTED TO GRAV. CALCULATION BETWEEN EACH OTHERS
'''ADDING TARGET OBJECT'S NOMINAL AND CLONES
ast_data = pd.read_csv(datadir)
for i in range(len(ast_data)):
    labl = 'A{:03d}'.format(ast_data['index'][i])
    asmd = ast_data['asmd'][i]
    eccs = ast_data['eccs'][i]
    incl = ast_data['incl'][i]
    ascn = ast_data['ascn'][i]
    peri = ast_data['peri'][i]
    sim.add(m=0., a=asmd, e=eccs, inc=incl, Omega=ascn, omega=peri, hash=labl)
print("Nominal orbit and its " + str(len(ast_data - 1)) + " clones of asteroid orbit has been added.")
sim.move_to_com()
'''

'''ADDING TARGET OBJECT'S NOMINAL ORBIT'''
ast_data = pd.read_csv(datadir + 'params_{a}.csv'.format(a=obj))
labl = 'A{:03d}'.format(ast_data['index'][0])
asmd = ast_data['asmd'][0]
eccs = ast_data['eccs'][0]
incl = ast_data['incl'][0]
ascn = ast_data['ascn'][0]
peri = ast_data['peri'][0]
sim.add(m=0., a=asmd, e=eccs, inc=incl, Omega=ascn, omega=peri, hash=labl)
print("Nominal orbit of asteroid has been added.")
sim.move_to_com()

'''SAVING INITIAL ORBITAL AND PHYSICAL PROPERTIES OF SOLAR SYSTEM OBJECTS AND NOMINAL ASTEROID IN CSV'''
pname = ['Sun', 'Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptunus', obj]
pmass = [ps[0].m, ps[1].m, ps[2].m, ps[3].m, ps[4].m, ps[5].m, ps[6].m, ps[7].m, ps[8].m, ps[9].m, ps['A000'].m]
peccn = [0, ps[1].e, ps[2].e, ps[3].e, ps[4].e, ps[5].e, ps[6].e, ps[7].e, ps[8].e, ps[9].e, ps['A000'].e]
pdist = [0, ps[1].a, ps[2].a, ps[3].a, ps[4].a, ps[5].a, ps[6].a, ps[7].a, ps[8].a, ps[9].a, ps['A000'].a]
pincl = [0, m.degrees(ps[1].inc), m.degrees(ps[2].inc), m.degrees(ps[3].inc), m.degrees(ps[4].inc), m.degrees(ps[5].inc), m.degrees(ps[6].inc), m.degrees(ps[7].inc), m.degrees(ps[8].inc), m.degrees(ps[9].inc), m.degrees(ps['A000'].inc)]
plasn = [0, m.degrees(ps[1].Omega), m.degrees(ps[2].Omega), m.degrees(ps[3].Omega), m.degrees(ps[4].Omega), m.degrees(ps[5].Omega), m.degrees(ps[6].Omega), m.degrees(ps[7].Omega), m.degrees(ps[8].Omega), m.degrees(ps[9].Omega), m.degrees(ps['A000'].Omega)]
pargp = [0, m.degrees(ps[1].omega), m.degrees(ps[2].omega), m.degrees(ps[3].omega), m.degrees(ps[4].omega), m.degrees(ps[5].omega), m.degrees(ps[6].omega), m.degrees(ps[7].omega), m.degrees(ps[8].omega), m.degrees(ps[9].omega), m.degrees(ps['A000'].omega)]
pmean = [0, m.degrees(ps[1].M), m.degrees(ps[2].M), m.degrees(ps[3].M), m.degrees(ps[4].M), m.degrees(ps[5].M), m.degrees(ps[6].M), m.degrees(ps[7].M), m.degrees(ps[8].M), m.degrees(ps[9].M), m.degrees(ps['A000'].M)]
dict = {"OBJ": pname, "MASS": pmass, "e": peccn, "a": pdist, "i": pincl, "Omega": plasn, "omega": pargp, "m": pmean}
di = pd.DataFrame(dict)
di.to_csv('grpn_' + start + '_init.csv')

'''INITIATE REBOUNDX AND ADDITIONAL FORCES'''
simx = reboundx.Extras(sim)
gr = simx.load_force(gr_type)
simx.add_force(gr)
gr.params["c"] = lightsp
ps[0].params["gr_source"] = 1

'''INITIATE AND WRITE THE SIMULATION LOGS'''
file = open("grpn_" + start + "_summary.txt","w")
file.write("SIMULATION SUMMARY LOGS \n")
file.write(start)
file.write("---------- \n")
file.write("\n")
file.write("TARGET OBJECT \n")
file.write("Object name                 : {a} \n".format(a=obj))
file.write("Epoch                       : {a} \n".format(a=date))
file.write("---------- \n")
file.write("\n")
file.write("SIMULATION SETTINGS \n")
file.write("Integrator                  : IAS15 \n")
file.write("Integrator tolerance        : {a} \n".format(a=sim.ri_ias15.epsilon))
file.write("Collision detection mode    : direct \n")
file.write("Collision resolver mode     : halt \n")
file.write("Simulation time limit (yr)  : {a} \n".format(a=sim_time))
file.write("Sim. data time step (yr)    : {a} \n".format(a=sim_dast))
file.write("---------- \n")

'''STARTING THE INTEGRATION AND DATA COLLECTION FOR EVERY TIMESTEP'''
print('NUMERICAL INTEGRATION IS IN PROGRESS')
n = 0
c = []                          # array for colliding clone registers
sim.heartbeat = heartbeat
sim.integrate(sim_time)  

'''ENDING THE SIMULATION PROCESS'''
print('-----')
print('SIMULATION IS FINISHED')
print("-----")

fn_time = datetime.now()
finish = fn_time.strftime("%Y%m%d_%H%M%S")
print(finish)

file.write("---------- \n")
file.write("The program started on {a} \n".format(a=start))
file.write("and ended on {a} \n".format(a=finish))
file.close()