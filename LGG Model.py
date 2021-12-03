# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

@author: Group 47
"""
import numpy as np

# COMBUSTION CHAMBER VALUES
D_c = 100e-3 # Diameter of the Combustion Chamber
L_c = 0.20   # Length of the Combustion Chamber
P_c0 = 101.3e+3  # Initial pressure in the Combustion Chamber (before detonation)
gamma_c = 1.4 # Gamma for the combustion products

# PISTON VALUES
m_pis = 600e-3 # Mass of the piston
mu_pis = 0.4   # Coefficient of friction for the piston against the pump tube

# PUMP TUBE VALUES
P_t0 = 300e+3  # Initial pressure in the pump tube (ahead of the piston)
L_t0 = 0.8     # Length of the pump tube
D_t = 12e-3   # Diameter of the pump tube

# RUPTURE DISK VALUES
P_rupt = 16e+6 # Pressure at which the rupture disk ruptures

# BARREL VALUES
L_b = 1   # Length of the barrel
D_b = 5e-3 # Diameter of the barrel
P_b0 = 0.1 # Initial upstream pressure in the barrel
gamma_b = 1.4 # Gamma for the light gas gun

# PROJECTILE & SABOT VALUES
m_p = 0.01e-3
m_s = m_p

# SIMULATION DETAILS
delta_t = 1e-5 # Time step length

# Calculating areas from diameters
A_c = np.pi * (D_c/2)**2
A_t = np.pi * (D_t/2)**2
A_b = np.pi * (D_b/2)**2


def DoIt(A_c = A_c, L_c = L_c, P_c0 = P_c0, gamma_c = gamma_c,\
         m_pis = m_pis, mu_pis = mu_pis, \
         P_t0 = P_t0, L_t0 = L_t0, A_t = A_t,\
         P_rupt = P_rupt,\
         L_b =  L_b, D_b = D_b, P_b0 = P_b0, gamma_b = gamma_b,\
         m_p = m_p, m_s = m_s,\
         delta_t = delta_t):
    # DEFINING ARRAYS
    P_t_array = []
    P_b_array = []
    x_pis_array = []
    x_p_array = []
    v_pis_array = []
    v_p_array = []
    a_pis_array = []
    a_p_array = []
    V_pis_array = []
    V_p_array = []
    
    
