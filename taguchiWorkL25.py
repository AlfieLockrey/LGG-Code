# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:04:35 2022

@author: euandh
"""

import lggmodel as lgg 
import numpy as np
from matplotlib import pyplot as plt

timeStep = 1e-6

def possibleCombinations(listOfLists):
    return [list(x) for x in np.array(np.meshgrid(*listOfLists)).T.reshape(-1,len(listOfLists))]

def areaCalc(d):
    return np.power(d/2, 2)*np.pi


# Constant values that were previously DPs
D_c = 19.685e-3
L0_c = 150e-3
D_pis = None  # equal to the D_pt
m_pr = 1.45e-3
C = 20e-3
D_b = 12.7e-3
gamma_lg = 1.41

def runFromArray(DP, NF):  
    """
    # This version uses the value of the piston diameter for the piston diameter, rather than basing it on the pt diameter
    return lgg.DoIt(D_c=DP[0], L0_c=DP[1], P0_c=NF[3], gamma_c=NF[0], D_pis=DP[2],
             mu_static_pis = 0 , mu_dynamic_pis = 0,
             P0_pt=DP[3], L0_pt=DP[6], D_pt=DP[7], P_rupt=DP[8], L_b=DP[9],
             D_b=DP[10], P0_b=NF[1], gamma_lg=DP[11], m_pr=DP[4], m_sb=0,
             mu_sb=0, gamma_ic=NF[2], V_ic=DP[5], C = DP[12], BR_exp = NF[4], 
             delta_t=timeStep, printing = False)
    """
    # This version matches the piston diameter to the pump tube diameter
    return lgg.DoIt(D_c=D_c, L0_c=L0_c, P0_c=NF[3], gamma_c=NF[0], D_pis=DP[3],
             mu_static_pis = 0 , mu_dynamic_pis = 0,
             P0_pt=DP[0], L0_pt=DP[2], D_pt=DP[3], P_rupt=DP[4], L_b=DP[5],
             D_b=D_b, P0_b=NF[1], gamma_lg=gamma_lg, m_pr=m_pr, m_sb=0,
             mu_sb=0, gamma_ic=NF[2], V_ic=DP[1], C = C, BR_exp = NF[4], 
             delta_t=timeStep, printing = False)

""" TAGUCHI ORTHAGONAL ARRAYS """
L12 = np.genfromtxt("taguchiL12.csv", delimiter = ',')
L12[0][0] = 1

L8 = np.genfromtxt("taguchiL8.csv", delimiter = ',')
L8[0][0] = 1

L25 = np.genfromtxt("taguchiL25.csv", delimiter = ',')
L25[0][0] = 1



""" DESIGN PARAMETERS  for L25"""
P0_pt = np.linspace(30e+6, 60e+6, 5)             #0
V_ic = [0.5, 0.75, 1, 1.25, 1.5]                   #1
L0_pt = [1, 1.25, 1.5, 1.75, 2]                    #2
D_pt = np.linspace(12.7e-3, 40e-3, 5)         #3
P_rupt = np.linspace(10e+6, 100e+6, 5)          #4
L_b = [1, 2, 3, 4, 5]                    #5
"""

      #     0   1       2      3     4      5     6     7     8       9   10
DP_list = [D_c, L0_c, D_pis, P0_pt, m_pr, V_ic, L0_pt, D_pt, P_rupt, L_b, D_b, gamma_lg, C]
DP_exp = possibleCombinations(DP_list)
"""
#           0       1    2       3       4     5
DP_list = [P0_pt, V_ic, L0_pt, D_pt, P_rupt, L_b]

""" NOISE FACTORS """
gamma_c = [1.2 , 1.3]     #0
P0_b = [100, 200]               #1
gamma_ic = [1.4, 1.4]           #2
P0_c = [101.3e+3, 101.3e+3]         #3
BR_exp = [0.80, 0.82]     #4

#             0       1       2       3      4
NF_list = [gamma_c, P0_b, gamma_ic, P0_c, BR_exp]
NF_exp = possibleCombinations(NF_list)

""" APPLY TAGUCHI'S METHOD """
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
expL25 = []
for i in range(0, len(L25)):
    expL25row = []
    for j in range(0, len(L25[i])):
        expL25row.append(DP_list[j][int(L25[i][j] - 1)])
    expL25.append(expL25row)
    
expL8 = []
for i in range(0, len(modL8)):
    expL8row = []
    for j in range(0, len(modL8[i])):
        expL8row.append(NF_list[j][int(modL8[i][j] - 1)])
    expL8.append(expL8row)    
        
# Run experiments
results = []
maxPresults = []
maxv_pisresults = []
def runExperiments():
    for i in range(0, len(expL25)):
        resultsRow = []
        maxPresultsrow = []
        maxv_pisresultsrow = []
        for j in range(0, len(expL8)):
            result = runFromArray(expL25[i], expL8[j])
            resultsRow.append(result[0])
            maxPresultsrow.append(result[1])
            maxv_pisresultsrow.append(result[2])
        results.append(resultsRow)
        maxPresults.append(maxPresultsrow)
        maxv_pisresults.append(maxv_pisresultsrow)
        print(i+1,"/",len(expL25))
    
    
""" SIGNAL TO NOISE RATIOS """
SNratios = []
def ratios():
    
    for i in range(0, len(results)):
        SNratio = 0
        for j in range(0, len(results[i])):
            SNratio += (1/len(DP_list))*(1/np.power(results[i][j],2))
        SNratios.append(-10*np.log10(SNratio))



DPs = ["P0_pt","V_ic","L0_pt","D_pt","P_rupt","L_b"]

def plotForDP(ind, plotAll = False):
    levels = DP_list[ind]
    averageRatios = []
    for level in levels:
        ratiosThisLevel = []
        for j in range(0, len(expL25)):
            if expL25[j][ind] == level and np.isnan(SNratios[j]) == False and\
                np.isinf(SNratios[j]) == False:
                ratiosThisLevel.append(SNratios[j])
            
        
        
        averageRatios.append(np.average(ratiosThisLevel))
    """
    print("---",DPs[i],"---")
    print("Max: ", max(averageRatios), " Min: ", min(averageRatios))
    print("Difference: ", max(averageRatios) - min(averageRatios))
    """
    diff = max(averageRatios) - min(averageRatios)  
    if plotAll is False:
        plt.plot(levels, averageRatios)
    if plotAll is True:
        plt.plot([1,2,3], averageRatios)
    return levels, averageRatios, diff


runExperiments()
ratios()

differences = []
avgRatios = []
for i in range(0, len(DP_list)):
    resultsOut = plotForDP(i, False)
    differences.append(resultsOut[2])
    avgRatios.append(resultsOut[1])
    plt.title(DPs[i])
    plt.xlabel("Levels")
    plt.ylabel("Average S/N Ratio")
    plt.figure()
"""
print("Max difference is ", max(differences))
"""
orderedSNratios = SNratios.copy()
orderedSNratios.sort(reverse=True)

"""
for i in range(0, len(orderedSNratios)):
    if np.isnan(orderedSNratios[i]) is True or np.isinf(orderedSNratios[i]) is True:
        orderedSNratios.remove(orderedSNratios[i])
"""

workable = False
index = -1

while workable is False:
    index+= 1
    L25index = SNratios.index(orderedSNratios[index])
    
    # Check if P0_pt < P_rupt and that the values are nan or +/- inf
    if expL25[L25index][0] >= expL25[L25index][4] or np.isnan(orderedSNratios[index]) == True or np.isinf(orderedSNratios[index]) == True:
        workable = False
        print(L25index, " is not workable")
    else:
        workable = True
        print(L25index, " is workable")

    bestDesign = SNratios.index(orderedSNratios[index])
    
print("\nBest design is experiment", bestDesign,"giving an average v =", np.average(results[bestDesign]))
print("\t Avg. Max P_pt [MPa] = ", np.average(maxPresults[bestDesign])/1e+6)
print("\t MPa, at most: [MPa]", max(maxPresults[bestDesign])/1e+6)
print("\t Avg. Max V_pis = ", np.average(maxv_pisresults[bestDesign]))
print("\t at most: ", max(maxv_pisresults[bestDesign]))
for i in range(0, len(expL25[bestDesign])):
    if i == 0 or i == 4:
        print(DPs[i]+":\t"+str(expL25[bestDesign][i]/1e+6) + " MPa")
    else:
        print(DPs[i]+":\t"+str(expL25[bestDesign][i]))
        
def exportResults():
    np.savetxt("expL8.csv", expL8, delimiter = ',')
    np.savetxt("expL25.csv", expL25, delimiter = ',')
    np.savetxt("resultsTaguchi.csv", results, delimiter = ',')
    
    L8_T = np.transpose(expL8)
    
    full = []
    for i in range(0, len(L8_T)):
        row = []
        for j in range(0, len(expL25[0])):
            row.append("")
        for j in range(0, len(L8_T[i])):
            row.append(L8_T[i][j])
        full.append(row)
    for i in range(0, len(expL25)):
        row = []
        for j in range(0, len(expL25[i])):
            row.append(expL25[i][j])
        for j in range(0, len(results[i])):
            row.append(results[i][j])
        row.append(SNratios[i])
        full.append(row)
    np.savetxt("TaguchiTable.csv", full, delimiter = ',',fmt = "%s")
    

    