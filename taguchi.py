# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:04:35 2022

@author: euandh
"""

import lggmodel as lgg 
import numpy as np

timeStep = 1e-6

def possibleCombinations(listOfLists):
    return [list(x) for x in np.array(np.meshgrid(*listOfLists)).T.reshape(-1,len(listOfLists))]

def areaCalc(d):
    return np.power(d/2, 2)*np.pi

def runFromArray(DP, NF):
    return lgg.DoIt(A_c = areaCalc(DP[0]), L0_c = DP[1], P0_c = NF[3], \
             gamma_c = NF[0], A_pis = areaCalc(DP[2]), mu_static_pis = 0,\
             mu_dynamic_pis = 0, P0_pt = DP[3], L0_pt = DP[6], \
             A_pt = areaCalc(DP[7]), P_rupt = DP[8], L_b = DP[9], D_b = DP[10],\
             P0_b = NF[1], gamma_lg = NF[4], m_pr = DP[4], m_sb = 0, \
             mu_sb = 0, gamma_ic = NF[2], V_ic = DP[5], delta_t = timeStep, \
             printing = False)

""" TAGUCHI ORTHAGONAL ARRAYS """
L12 = np.genfromtxt("taguchiL12.csv", delimiter = ',')
L12[0][0] = 1

L8 = np.genfromtxt("taguchiL8.csv", delimiter = ',')
L8[0][0] = 1

""" DESIGN PARAMETERS """
D_c = [0.01, 0.05]
L0_c = [5e-3, 30e-3]
D_pis = [5e-3, 5e-2]
P0_pt = [5e+6, 5e+7]
m_pr = [1e-3, 5e-2]
V_ic = [1, 1.5]
L0_pt = [0.75, 1.25]
D_pt = [5e-3, 1e-1]
P_rupt = [1e+6, 5e+7]
L_b = [0.75, 1.25]
D_b = [1e-3, 5e-2]

      #     0   1       2      3     4      5     6     7     8       9   10
DP_list = [D_c, L0_c, D_pis, P0_pt, m_pr, V_ic, L0_pt, D_pt, P_rupt, L_b, D_b]
DP_exp = possibleCombinations(DP_list)

""" NOISE FACTORS """
gamma_c = [1.2, 1.4]
P0_b = [100, 200]
gamma_ic = [1.35, 1.45]
P0_c = [90e+3, 110e+3]
gamma_lg = [1.41, 1.667]

#             0       1       2       3      4
NF_list = [gamma_c, P0_b, gamma_ic, P0_c, gamma_lg]
NF_exp = possibleCombinations(NF_list)

"""APPLY TAGUCHI"""
# Modify L8
modL8 = []
for row in L8:
    modRow = []
    for element in row:
        if element == 3:
            modRow.append(1.0)
        elif element == 4:
            modRow.append(2.0)
        else:
            modRow.append(element)
    modL8.append(modRow)
    
# Use orthogonal matrices to set levels for experiments
expL12 = []
for i in range(0, len(L12)):
    expL12row = []
    for j in range(0, len(L12[i])):
        expL12row.append(DP_list[j][int(L12[i][j] - 1)])
    expL12.append(expL12row)
    
expL8 = []
for i in range(0, len(modL8)):
    expL8row = []
    for j in range(0, len(modL8[i])):
        expL8row.append(NF_list[j][int(modL8[i][j] - 1)])
    expL8.append(expL8row)    
        
# Run experiments
results = []
for i in range(0, len(L12)):
    resultsRow = []
    for j in range(0, len(L8)):
        resultsRow.append(runFromArray(L12[i], L8[j]))
    results.append(resultsRow)
    print(i+1,"/",len(L12))
