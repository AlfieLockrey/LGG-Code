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

def runFromArray(DP, NF):  
    return lgg.DoIt(D_c=DP[0], L0_c=DP[1], P0_c=NF[3], gamma_c=NF[0], D_pis=DP[2],
             mu_static_pis = 0 , mu_dynamic_pis = 0,
             P0_pt=DP[3], L0_pt=DP[6], D_pt=DP[7], P_rupt=DP[8], L_b=DP[9],
             D_b=DP[10], P0_b=NF[1], gamma_lg=DP[11], m_pr=DP[4], m_sb=0,
             mu_sb=0, gamma_ic=NF[2], V_ic=DP[5], C = DP[12], BR_exp = NF[4], 
             delta_t=timeStep, printing = False)

""" TAGUCHI ORTHAGONAL ARRAYS """
L12 = np.genfromtxt("taguchiL12.csv", delimiter = ',')
L12[0][0] = 1

L8 = np.genfromtxt("taguchiL8.csv", delimiter = ',')
L8[0][0] = 1

L27 = np.genfromtxt("taguchiL27.csv", delimiter = ',')
L27[0][0] = 1

""" DESIGN PARAMETERS """
D_c = [0.025, 0.03, 0.045]             #0
L0_c = [150e-3, 200e-3, 250e-3]        #1
D_pis = [5e-3, 30e-3, 70e-3]           #2
P0_pt = [1e+5, 1e+6, 1e+7]             #3
m_pr = [5e-3, 15e-3, 30e-3]            #4
V_ic = [0.5, 1, 1.5]                   #5
L0_pt = [1, 1.5, 2]                    #6
D_pt = [30e-3, 45e-3, 60e-3]           #7
P_rupt = [1e+6, 5e+6, 10e+6]           #8
L_b = [1.25, 1.75, 2.25]               #9
D_b = [6.35e-3, 12.7e-3, 19.05e-3]     #10
gamma_lg = [1.41, 1.667, 1.667]        #11
C = [15e-3, 20e-3, 25e-3]              #12

      #     0   1       2      3     4      5     6     7     8       9   10
DP_list = [D_c, L0_c, D_pis, P0_pt, m_pr, V_ic, L0_pt, D_pt, P_rupt, L_b, D_b, gamma_lg, C]
DP_exp = possibleCombinations(DP_list)

""" NOISE FACTORS """
gamma_c = [1.2 , 1.4]     #0
P0_b = [100, 200]               #1
gamma_ic = [1.2, 1.4]           #2
P0_c = [90e+3, 110e+3]         #3
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
expL27 = []
for i in range(0, len(L27)):
    expL27row = []
    for j in range(0, len(L27[i])):
        expL27row.append(DP_list[j][int(L27[i][j] - 1)])
    expL27.append(expL27row)
    
expL8 = []
for i in range(0, len(modL8)):
    expL8row = []
    for j in range(0, len(modL8[i])):
        expL8row.append(NF_list[j][int(modL8[i][j] - 1)])
    expL8.append(expL8row)    
        
# Run experiments
results = []
def runExperiments():
    
    dicts = []
    for i in range(0, len(expL27)):
        resultsRow = []
        for j in range(0, len(expL8)):
            resultsRow.append(runFromArray(expL27[i], expL8[j]))
            dicts.append({
                         
                         })
        results.append(resultsRow)
        print(i+1,"/",len(expL27))
    
    
""" SIGNAL TO NOISE RATIOS """
SNratios = []
def ratios():
    
    for i in range(0, len(results)):
        SNratio = 0
        for j in range(0, len(results[i])):
            SNratio += (1/len(DP_list))*(1/np.power(results[i][j],2))
        SNratios.append(-10*np.log10(SNratio))

def exportResults():
    np.savetxt("expL8.csv", expL8, delimiter = ',')
    np.savetxt("expL27.csv", expL27, delimiter = ',')
    np.savetxt("resultsTaguchi.csv", results, delimiter = ',')

DPs = ["D_c","L0_c","D_pis","P0_pt","m_proj","V_ic","L0_pt","D_pt","P_rupt","L_b","D_b","gamma_lg","C"]

def plotForDP(ind, plotAll = False):
    levels = DP_list[ind]
    averageRatios = []
    for level in levels:
        ratiosThisLevel = []
        for j in range(0, len(expL27)):
            if expL27[j][ind] == level and np.isnan(SNratios[j]) == False and\
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
    
    if P0_pt >= P_rupt or np.isnan(orderedSNratios[index]) == True or np.isinf(orderedSNratios[index]) == True:
        workable = False
    else:
        workable = True

    bestDesign = SNratios.index(orderedSNratios[index])
    
print("\nBest design is experiment", bestDesign,"giving an average v =", np.average(results[bestDesign]))
for i in range(0, len(expL27[bestDesign])):
    print(DPs[i], ":\t", expL27[bestDesign][i])
    