# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:54:23 2022

@author: euandh
"""

import numpy as np
import lggmodel as lgg

""" 0. RUN VALUES """
numberOfPoints = 1000
samples = 30
NFNoOfPoints = numberOfPoints
deltaT = 1e-6
secondsPerPointEst = 0.034  #Estimate of how long it takes to do a single run

""" ------------ 1. VARIABLES ------------ """
""" 1A - Constant Values """
D_c = 19.685e-3
L0_c = 150e-3
# D_pis = None  # equal to the D_pt
m_pr = 1.45e-3
C = 20e-3
D_b = 12.7e-3
gamma_lg = 1.41
gamma_ic = 1.4
P0_c = 101.3e+3

""" 1B - Design Parameter Ranges """
P0_pt = np.linspace(30e+6, 60e+6, numberOfPoints)       #0
V_ic = np.linspace(0.5, 1.5, numberOfPoints)            #1
L0_pt = np.linspace(1, 2, numberOfPoints)               #2
D_pt = np.linspace(12.7e-3, 40e-3, numberOfPoints)      #3 
P_rupt = np.linspace(10e+6, 100e+6, numberOfPoints)     #4 
L_b = np.linspace(1, 5, numberOfPoints)                 #5

""" 1C - Noise Factor Ranges """
gamma_c = 1.3
P0_b = 101.3e+3
BR_exp = 0.818

""" ------------ 2. ALL COMBINATIONS (ITERATION) ------------ """
def numberOfCombinations(listOfLists):
    # Returns the number of combinations that there are
    m = len(listOfLists)
    s = len(listOfLists[0])
    return s**m

def possibleCombinations(listOfLists):
    # Returns an array with all the possible different combinations
    return [list(x) for x in np.array(np.meshgrid(*listOfLists)).T.reshape(-1,len(listOfLists))]

def areaCalc(d):
    return np.power(d/2, 2)*np.pi

listOfLists = [P0_pt, V_ic, L0_pt, D_pt, P_rupt, L_b]
combCount = numberOfCombinations(listOfLists)
print("Number of Combinations: "+str(combCount))
print("Estimated Time: "+str((combCount*secondsPerPointEst)/(60*60))+" hours")

def runAllCombinations(listOfLists):
    combCount = numberOfCombinations(listOfLists)
    print("Number of combinations = "+str(combCount))
    combinations = possibleCombinations(listOfLists)
    for i in range(0, len(combinations)):
        print(str(i + 1)+"/"+str(combCount))
        results = lgg.DoIt(D_c=D_c, L0_c=L0_c, P0_c=P0_c, gamma_c=gamma_c, D_pis=combinations[i][3],
                 mu_static_pis=0, mu_dynamic_pis=0,
                 P0_pt=combinations[i][0], L0_pt=combinations[i][2], D_pt=combinations[i][3], P_rupt=combinations[i][4], L_b=combinations[i][5],
                 D_b=D_b, P0_b=P0_b, gamma_lg=gamma_lg, m_pr=m_pr, m_sb=0,
                 mu_sb=0, gamma_ic=gamma_ic, V_ic=combinations[i][1], delta_t=deltaT)
        for j in range(0, len(results)):
            combinations[i].append(results[j])
        
    return combinations

        

