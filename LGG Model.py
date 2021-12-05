# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 11:44:15 2021

@author: Group 47
"""
import numpy as np
from matplotlib import pyplot as plt

# TROUBLESHOOTING
Steps = 1000

# COMBUSTION CHAMBER VALUES
D_c = 100e-3 # Diameter of the Combustion Chamber
L_c = 0.20   # Length of the Combustion Chamber
P_c0 = 5e+6  # Initial pressure in the Combustion Chamber (before detonation)
gamma_c = 1.4 # Gamma for the combustion products

# PISTON VALUES
D_pis = 12e-3 
m_pis = 600e-3 # Mass of the piston
mu_pis = 0.4   # Coefficient of friction for the piston against the pump tube

# PUMP TUBE VALUES
P_t0 = 300e+3  # Initial pressure in the pump tube (ahead of the piston)
L_t0 = 0.8     # Length of the pump tube
D_t = D_pis  # Diameter of the pump tube

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
delta_t = 1e-4 # Time step length

# Calculating areas from diameters
A_c = np.pi * (D_c/2)**2
A_t = np.pi * (D_t/2)**2
A_b = np.pi * (D_b/2)**2
A_pis = np.pi * (D_pis/2)**2


def DoIt(A_c = A_c, L_c = L_c, P_c0 = P_c0, gamma_c = gamma_c,\
         A_pis = A_pis, m_pis = m_pis, mu_pis = mu_pis, \
         P_t0 = P_t0, L_t0 = L_t0, A_t = A_t,\
         P_rupt = P_rupt,\
         L_b =  L_b, D_b = D_b, P_b0 = P_b0, gamma_b = gamma_b,\
         m_p = m_p, m_s = m_s,\
         delta_t = delta_t):
    # DEFINING ARRAYS
    """
    any initial value that's in one of these arrays that is NOT a 
    variable will need to be replaced, I've just put in some placeholder
    values so that I can test the program.
    - Euan
    """
    P_c_array = [P_c0]
    P_t_array = [P_t0]
    P_b_array = []
    x_pis_array = [0.01]
    x_p_array = [0]
    v_pis_array = [0]
    v_p_array = [0]
    a_pis_array = [0]
    a_p_array = [0]
    V_pis_array = [0]
    V_p_array = [0]
    n_array = [] # Counts the time step number
    t_array = [] # Holds the time step value
    
    # COMBUSTION CHAMBER
        # come back to later lol
    
    diskBroken = False
    
    def pistonElement():
        delta_P = P_c_array[-1] - P_t_array[-1]
        F_pressure = delta_P*A_pis
        F_fric = mu_pis*0  # change later
        F_res = F_pressure - F_fric
        a_pis = F_res/m_pis
        a_pis_array.append(a_pis)
        v_pis = v_pis_array[-1] + a_pis*delta_t
        x_pis = v_pis_array[-1]*delta_t + 0.5*a_pis*(delta_t**2) + x_pis_array[-1] 
        v_pis_array.append(v_pis)
        x_pis_array.append(x_pis)
        
        
        # Not entirely sure this is actually how we should do this:
        P_c = P_c_array[0]*np.power(x_pis_array[0]/x_pis_array[-1], gamma_c)
        P_c_array.append(P_c)
        # Ahead pressure, this line is causing issues
        P_t = P_t_array[0]*np.power((L_t0 - x_pis_array[0])/(L_t0 - x_pis_array[-1]), gamma_b)
        P_t_array.append(P_t)
    
    i = 0
    n_array.append(i)
    t_array.append(i*delta_t)
    #while x_p_array[-1] < L_b:
    while n_array[-1] < Steps:
        if P_t_array[-1] < P_rupt and diskBroken == False:
            pistonElement()
            i += 1
            n_array.append(i)
            t_array.append(i*delta_t)
            print(P_t_array[-1])
        else:
            print("Barrel Reached")
            x_p_array.append(L_b)
    
    return P_t_array, P_c_array, x_pis_array, a_pis_array, n_array, t_array

Steps = 1e+4
data = DoIt()
fig, ax = plt.subplots(constrained_layout=True)
ax.plot(data[4], data[2], color = "blue")
ax.tick_params(axis = 'y', labelcolor = "blue")
ax.set_xlabel("Step number")
ax.set_ylabel("Position [m]", color = "blue")
ax2 = ax.twinx()
ax2.plot(data[4], data[3], color =  "red")
ax2.tick_params(axis = 'y', labelcolor = "red")
ax2.set_ylabel("Acceleration [m/s^2]", color = "red")
plt.show()