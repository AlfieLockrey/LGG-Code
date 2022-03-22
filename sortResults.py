# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 11:36:33 2022

@author: euandh
"""

""" CAPPED VALUES """
P_pt_max = 250e+6
P_cc_max = 300e+6
# Piston_x_max can't be beyond the pump tube length (whoops I didn't account for this in the code)
V_pis = 500
P_rupt = 150e+6 # HAS TO BE LARGER THAN THIS


results = [] # Imported from iterative or MC code
""" DATA PROCESSING """

def meetsConditions(result):
    

for i in range(0, len(results)):
    print("d")
